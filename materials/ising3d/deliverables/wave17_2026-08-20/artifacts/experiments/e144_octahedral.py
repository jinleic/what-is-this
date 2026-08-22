#!/usr/bin/env python3
"""Octahedral (order-48 point-group) feasibility-then-execution for C(3,2,6).

Wave 13/14 certified the extensive-charge class quotient

    Q(C) = rank S_C - rank M_C - 1 = 1

for the Z^3 radius-2 support-size classes C(3,2,s), s <= 4 (wave 13) and
s <= 5 (wave 14, `experiments/e135_extensive_3d_s5.py`), by an exact
column echelon over F_p at two primes against the analytic {I, h} kernel
bound rank_Q M_C <= rank S_C - 2.  The s <= 6 class was NOT attempted.
This script answers the wave-15 question: can the octahedral point group
(the 48 signed coordinate permutations, the point group of the Z^3 lattice
preserving the 3x3x3 box about its centre) compress the 21,121,156-column
s <= 5 system -- and the 236,912,446-column s <= 6 system -- to orbit
representatives, and can that compression DECIDE the s <= 6 quotient?

Method (feasibility first, all exact integer arithmetic):

* FEASIBILITY.  |C(3,2,s)| = sum_{k<=s} C(27,k) 3^k by the closed form.
  The number of G-orbits on the word set is computed by the Burnside
  lemma: for each of the 48 site permutations the fixed words are exactly
  the letter assignments constant on the site cycles, counted by a DP over
  the cycle multiset (exact integers).  The SAME is done for the ROW side
  (anchored words = translation-orbit representatives, the rank-S rows):
  the induced action fixes an anchored word w iff g.w = tau_t . w for a
  unique t, i.e. iff the support is invariant under the affine map
  A_t(y) = g y - t with letters constant on A_t-orbits; counted by a DP
  over the finite in-box cycles of A_t (t ranges over a 125-point box).
  Validations: (a) at (g = e, t = 0) the DP reproduces the rank-S closed
  form 1 + sum_k A_k 3^k for every s; (b) independent brute-force orbit
  counting by canonical forms agrees at s <= 3 on the 27-site box for
  both families; (c) the 8-site box B_1 (half-integer centre) is checked
  by full brute force, columns through s = 8 and rows through s = 4.

* EXECUTION.  A G-orbit-representative column echelon (canonical word =
  orbit minimum under the lexicographic (support mask, letters) order)
  over the same two primes as wave 14.
  [THEOREM-elementary CAP] an echelon fed only orbit representatives
  certifies at most #orbits pivots,
  and the s <= 6 certificate needs rank S_6 - 2 = 191,949,608 pivots, so
  the representative echelon can only ever certify 5,048,368/191,949,608
  = 2.63% of the target: the s <= 6 quotient CANNOT be decided by column
  orbit compression.  The representative echelon is nevertheless run at
  s <= 5 (where the exact full answer 14,757,410 is known: the empirical
  stall fraction) and at s <= 6 (a genuine certified lower bound
  rank_Q M_6 >= rank_Fp(M_6 restricted to representatives)); the full-set
  certificate itself is walled out by exact memory models (lead table +
  CSR pivot store, both calibrated on the frozen wave-14 measurements).

Reproduce:

PYTHONPATH=src .venv/bin/python experiments/e144_octahedral.py
PYTHONPATH=src .venv/bin/python tests/test_octahedral.py
"""

from __future__ import annotations

import datetime as _dt
import gc
import hashlib
import itertools
import json
import math
import platform
import resource
import sys
import time
from array import array
from pathlib import Path
from typing import Iterator

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "integrability" / "octahedral_s6.json"
FROZEN_S5 = ROOT / "results" / "integrability" / "extensive_3d_s5.json"
SCRIPT = "experiments/e144_octahedral.py"

PRIMES = (2_147_483_647, 2_147_483_629)

# Predeclared budgets: process_time seconds (never wall clock), RSS bytes.
RSS_WALL_BYTES = 8 * 1024**3  # the wave-15 decision gate: 8 GiB
PREFLIGHT_BUDGET_S = 1200.0
S2_PRIME_BUDGET_S = 120.0
S3_PRIME_BUDGET_S = 900.0
S4_PRIME_BUDGET_S = 2400.0
REP5_PRIME_BUDGET_S = 1200.0
REP6_PRIME_BUDGET_S = 7200.0

# Frozen wave-13/14 certified anchors (results/integrability/extensive_3d_*.json).
CERTIFIED_S4_COLUMNS = 1_503_766
CERTIFIED_S4_RANK_S = 822_334
CERTIFIED_S4_RANK_M = 822_332
CERTIFIED_S5_COLUMNS = 21_121_156
CERTIFIED_S5_RANK_S = 14_757_412
CERTIFIED_S5_RANK_M = 14_757_410
# Frozen wave-14 headline measurements used by the memory models.
FROZEN_S5_PIVOTS = 14_757_410
FROZEN_S5_PIVOT_NNZ = 313_579_672  # excluding lead entries
FROZEN_S5_PEAK_RSS = 2_404_089_856
FROZEN_S5_SECONDS = 1_446.725794

GROUP_ORDER = 48

# Lead-table capacities (slots) and pivot guards per stage; guards leave
# substantial unused capacity below the table-full failure mode.
STAGE_PLAN = {
    # stage: (cap_bits, guard_pivots, budget_s per prime, log_interval)
    "full_s2": (11, 700, S2_PRIME_BUDGET_S, 1_000),
    "full_s3": (16, 33_000, S3_PRIME_BUDGET_S, 20_000),
    "full_s4": (21, 905_000, S4_PRIME_BUDGET_S, 250_000),
    "rep_s5": (20, 480_000, REP5_PRIME_BUDGET_S, 100_000),
    "rep_s6": (23, 5_100_000, REP6_PRIME_BUDGET_S, 500_000),
}


def peak_rss_bytes() -> int:
    """Darwin reports ru_maxrss in bytes; Linux reports KiB."""
    raw = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return raw if platform.system() == "Darwin" else raw * 1024


def process_time() -> float:
    """The only clock this script gates on (never wall clock, no alarms)."""
    return time.process_time()


class ResourceWall(RuntimeError):
    """An OBSERVED budget breach; carries the exact measured record."""

    def __init__(self, stage: str, reason: str, processed: int, total: int,
                 elapsed: float, rss: int) -> None:
        super().__init__(reason)
        self.record = {
            "claim_tag": "[COMPUTATION]",
            "kind": "observed_resource_wall",
            "stage": stage,
            "reason": reason,
            "processed": processed,
            "total": total,
            "elapsed_seconds_process_time": round(elapsed, 6),
            "peak_rss_bytes": rss,
            "peak_rss_mib": round(rss / 1024**2, 3),
        }


# ---------------------------------------------------------------------------
# Finite geometry: the 3x3x3 box B_2 inside the extended box C_2 (wave-14
# convention, identical to experiments/e135_extensive_3d_s5.py).
# ---------------------------------------------------------------------------

BASE_SITES: tuple[tuple[int, int, int], ...] = tuple(itertools.product(range(3), repeat=3))
BASE_INDEX = {site: i for i, site in enumerate(BASE_SITES)}
EXT_SITES: tuple[tuple[int, int, int], ...] = tuple(itertools.product(range(-1, 4), repeat=3))
EXT_INDEX = {site: i for i, site in enumerate(EXT_SITES)}
BASE_BITS = tuple(1 << EXT_INDEX[site] for site in BASE_SITES)
EXT_MASK = (1 << len(EXT_SITES)) - 1
NEIGHBORS: dict[int, tuple[int, ...]] = {}
for _site, _bit in zip(BASE_SITES, BASE_BITS):
    _nearby = []
    for _axis in range(3):
        for _sign in (-1, 1):
            _nb = list(_site)
            _nb[_axis] += _sign
            _nearby.append(1 << EXT_INDEX[tuple(_nb)])
    NEIGHBORS[_bit] = tuple(_nearby)


def anchor_key(x_mask: int, z_mask: int) -> int:
    """Injective int64-safe key of the translation orbit of X^x Z^z.

    Wave-14 convention (tightened radix for the s <= 6 classes): the
    canonical representative is the unique translate whose support has
    coordinatewise minimum 0; relative coordinates dx, dy, dz lie in 0..4
    and the letter code in 1..3, so the per-site digit
    ((dx*5+dy)*5+dz)*3 + (code-1) lies in 0..374 < 375.  A word of
    support size n has n sorted digits, so packed < 375^7 < 2^60 for the
    largest reachable support (smax+1 = 7), and key = (size << 60) |
    packed is injective, fits a signed int64 (max < 7*2^60 + 375^7 <
    2^63), and orders lexicographically as (size, packed digits).
    """
    support = x_mask | z_mask
    if not support:
        return 0
    mnx = mny = mnz = 99
    coords = []
    m = support
    size = 0
    while m:
        bit = m & -m
        m ^= bit
        pos = bit.bit_length() - 1
        sx, sy, sz = EXT_SITES[pos]
        if sx < mnx:
            mnx = sx
        if sy < mny:
            mny = sy
        if sz < mnz:
            mnz = sz
        coords.append((pos, sx, sy, sz))
        size += 1
    vals = []
    for pos, sx, sy, sz in coords:
        code = (1 if (x_mask >> pos) & 1 else 0) | (2 if (z_mask >> pos) & 1 else 0)
        vals.append((((sx - mnx) * 5 + (sy - mny)) * 5 + (sz - mnz)) * 3 + (code - 1))
    vals.sort()
    packed = 0
    for v in vals:
        packed = packed * 375 + v
    return (size << 60) | packed


def image_terms(x_mask: int, z_mask: int) -> dict[tuple[int, int], int]:
    """Raw image of X^x Z^z under D = (1/2)[H, .], a=b=1, ordered X^a Z^b.

    Field:  (1/2)[X_s, X^a Z^b] = [b_s=1] X^(a+e_s) Z^b              (+1 each)
    Bond:   (1/2)[Z_uZ_v, X^a Z^b] = -[a_u xor a_v=1] X^a Z^(b+e_u+e_v)
    with u ranging over X-sites and v over non-X neighbours of u.
    """
    out: dict[tuple[int, int], int] = {}
    m = z_mask
    while m:
        bit = m & -m
        m ^= bit
        k = (x_mask ^ bit, z_mask)
        out[k] = out.get(k, 0) + 1
    m = x_mask
    while m:
        bit = m & -m
        m ^= bit
        for nb in NEIGHBORS[bit]:
            if not x_mask & nb:
                k = (x_mask, z_mask ^ bit ^ nb)
                out[k] = out.get(k, 0) - 1
    return {k: v for k, v in out.items() if v}


def class_column_count(smax: int) -> int:
    """Closed form: sum_{k<=smax} C(27,k) 3^k (each site carries X, Z or XZ)."""
    return sum(math.comb(27, k) * 3**k for k in range(smax + 1))


def iter_class_columns(smax: int) -> Iterator[tuple[int, int]]:
    """Stream the raw class columns: identity first, then by support size."""
    yield (0, 0)
    for size in range(1, smax + 1):
        for positions in itertools.combinations(range(27), size):
            bits = [BASE_BITS[p] for p in positions]
            for codes in itertools.product((1, 2, 3), repeat=size):
                x = z = 0
                for b, code in zip(bits, codes):
                    if code & 1:
                        x |= b
                    if code & 2:
                        z |= b
                yield (x, z)


def anchored_subset_count_inclusion_exclusion(k: int) -> int:
    """A_k = C(27,k) - 3 C(18,k) + 3 C(12,k) - C(8,k)."""
    return (math.comb(27, k) - 3 * math.comb(18, k)
            + 3 * math.comb(12, k) - math.comb(8, k))


def anchored_subset_count_bruteforce(k: int) -> int:
    n = 0
    for combo in itertools.combinations(range(27), k):
        sites = [BASE_SITES[c] for c in combo]
        if (min(s[0] for s in sites) == 0
                and min(s[1] for s in sites) == 0
                and min(s[2] for s in sites) == 0):
            n += 1
    return n


_ANCHORED_SUBSETS_CACHE: dict[int, int] = {}


def anchored_subset_count(k: int) -> int:
    """A_k with the inclusion-exclusion and brute-force values forced equal."""
    if k not in _ANCHORED_SUBSETS_CACHE:
        ie = anchored_subset_count_inclusion_exclusion(k)
        assert ie == anchored_subset_count_bruteforce(k), (k, ie)
        _ANCHORED_SUBSETS_CACHE[k] = ie
    return _ANCHORED_SUBSETS_CACHE[k]


def rank_S_closed_form(smax: int) -> int:
    """rank S_C = 1 + sum_k A_k 3^k (distinct anchored orbits of class words)."""
    return 1 + sum(anchored_subset_count(k) * 3**k for k in range(1, smax + 1))


# ---------------------------------------------------------------------------
# The octahedral group: 48 signed coordinate permutations about the box
# centre, acting on the 27 sites.  Doubled centred coordinates keep every
# action integral for both odd and even box sides.
# ---------------------------------------------------------------------------

GROUP_SIGNS: tuple[tuple[tuple[int, ...], tuple[int, ...]], ...] = tuple(
    (p, s)
    for p in itertools.permutations(range(3))
    for s in itertools.product((1, -1), repeat=3)
)
CTR = tuple((s[0] - 1, s[1] - 1, s[2] - 1) for s in BASE_SITES)  # centred coords


def _perm_table27(p: tuple[int, ...], signs: tuple[int, ...]) -> tuple[int, ...]:
    table = []
    for idx in range(27):
        c = CTR[idx]
        img = tuple(signs[i] * c[p[i]] + 1 for i in range(3))
        table.append(BASE_INDEX[img])
    return tuple(table)


PERMS27 = tuple(_perm_table27(p, s) for p, s in GROUP_SIGNS)


def cycle_lengths(perm: tuple[int, ...]) -> list[int]:
    seen = [False] * 27
    out = []
    for i in range(27):
        if seen[i]:
            continue
        j, length = i, 0
        while not seen[j]:
            seen[j] = True
            j = perm[j]
            length += 1
        out.append(length)
    return out


def fix_counts_from_cycles(cycles: list[int], smax: int) -> list[int]:
    """fix[s] = # words with support <= s invariant under the site permutation.

    A word is fixed iff every site cycle carries one letter (I or one of
    {X, Z, XZ}); a taken cycle of length L contributes L to the support
    size and a factor 3.
    """
    ways = [1] + [0] * smax
    for length in cycles:
        new = ways[:]
        for n in range(length, smax + 1):
            new[n] += 3 * ways[n - length]
        ways = new
    fix: list[int] = []
    total = 0
    for n in range(smax + 1):
        total += ways[n]
        fix.append(total)
    return fix


def column_orbit_counts(smax: int) -> tuple[list[int], list[list[int]]]:
    """Burnside: G-orbits on the word class C(3,2,s) for s = 0..smax."""
    fixes = [fix_counts_from_cycles(cycle_lengths(p), smax) for p in PERMS27]
    orbits: list[int] = []
    for s in range(smax + 1):
        total = sum(f[s] for f in fixes)
        if total % GROUP_ORDER:
            raise AssertionError(f"column Burnside not integral at s={s}: {total}")
        orbits.append(total // GROUP_ORDER)
    return orbits, fixes


# --- row side: induced action on anchored words, affine A_t(y) = g y - t ----

def affine_in_box_cycles(p: tuple[int, ...], signs: tuple[int, ...],
                         t: tuple[int, int, int]) -> list[frozenset[int]]:
    """Finite cycles of A_t (on centred coords) fully inside the 27-site box.

    A support set S is A_t-invariant iff S is a union of full A_t cycles;
    cycles that leave the box can never be part of an in-box invariant
    support, so only the in-box finite cycles matter.
    """
    cycles: list[frozenset[int]] = []
    seen_cycles: set[frozenset[int]] = set()
    done: set[int] = set()
    for start in range(27):
        if start in done:
            continue
        walk = [start]
        posmap = {start: 0}
        cur = start
        while True:
            cy = CTR[cur]
            ny = tuple(signs[i] * cy[p[i]] - t[i] for i in range(3))
            if not all(-1 <= v <= 1 for v in ny):
                break  # orbit leaves the box: start is not on an in-box cycle
            nxt = BASE_INDEX[tuple(v + 1 for v in ny)]
            if nxt in posmap:
                cyc = frozenset(walk[posmap[nxt]:])
                if cyc not in seen_cycles:
                    seen_cycles.add(cyc)
                    cycles.append(cyc)
                break
            posmap[nxt] = len(walk)
            walk.append(nxt)
            cur = nxt
        done.update(walk)
    return cycles


def row_fix_from_cycles(cycles: list[frozenset[int]], smax: int) -> list[int]:
    """Anchored words invariant under A_t with letter-constant cycles.

    dp[n][cov] counts letter assignments of exact support size n whose
    support union touches the face set cov (bit i = centred coordinate
    -1 present on axis i, i.e. the anchored min-face condition).  Each
    taken cycle carries one of 3 letters.  Returns cumulative counts for
    support sizes 1..s (the empty word is handled by the caller).
    """
    dp = [[0] * 8 for _ in range(smax + 1)]
    dp[0][0] = 1
    for cyc in cycles:
        length = len(cyc)
        cov = 0
        for i in cyc:
            c = CTR[i]
            if c[0] == -1:
                cov |= 1
            if c[1] == -1:
                cov |= 2
            if c[2] == -1:
                cov |= 4
        new = [row[:] for row in dp]
        for n in range(length, smax + 1):
            base = dp[n - length]
            for c0 in range(8):
                if base[c0]:
                    new[n][c0 | cov] += 3 * base[c0]
        dp = new
    return [sum(dp[n][7] for n in range(1, s + 1)) for s in range(smax + 1)]


def row_orbit_counts(smax: int) -> tuple[list[int], list[list[int]]]:
    """Burnside: G-orbits on anchored class words (the rank-S rows)."""
    per_g: list[list[int]] = []
    for p, signs in GROUP_SIGNS:
        fix = [1] * (smax + 1)  # the empty word is fixed by every element
        for t in itertools.product(range(-2, 3), repeat=3):
            dp = row_fix_from_cycles(affine_in_box_cycles(p, signs, t), smax)
            for s in range(smax + 1):
                fix[s] += dp[s]
        per_g.append(fix)
    orbits: list[int] = []
    for s in range(smax + 1):
        total = sum(f[s] for f in per_g)
        if total % GROUP_ORDER:
            raise AssertionError(f"row Burnside not integral at s={s}: {total}")
        orbits.append(total // GROUP_ORDER)
    return orbits, per_g


# ---------------------------------------------------------------------------
# Independent brute-force orbit counting (validation; canonical forms).
# ---------------------------------------------------------------------------

def brute_column_orbits(smax: int) -> int:
    canon: set[tuple[tuple[int, int], ...]] = set()
    for k in range(smax + 1):
        for sites in itertools.combinations(range(27), k):
            for codes in itertools.product((1, 2, 3), repeat=k):
                w = tuple(zip(sites, codes))
                best = None
                for perm in PERMS27:
                    key = tuple(sorted((perm[i], c) for i, c in w))
                    if best is None or key < best:
                        best = key
                canon.add(best)
    return len(canon)


def _reanchored_key(pairs: list[tuple[int, int]]) -> tuple:
    pts = [(BASE_SITES[i], c) for i, c in pairs]
    shift = [min(p[0][a] for p in pts) for a in range(3)]
    return tuple(sorted(
        ((p[0][0] - shift[0], p[0][1] - shift[1], p[0][2] - shift[2]), p[1])
        for p in pts
    ))


def brute_row_orbits(smax: int) -> int:
    canon: set[tuple] = set()
    for k in range(0, smax + 1):
        for sites in itertools.combinations(range(27), k):
            if k:
                pts = [BASE_SITES[i] for i in sites]
                if not (min(p[0] for p in pts) == 0
                        and min(p[1] for p in pts) == 0
                        and min(p[2] for p in pts) == 0):
                    continue
            for codes in itertools.product((1, 2, 3), repeat=k):
                w = tuple(zip(sites, codes))
                best = None
                for perm in PERMS27:
                    g = [(perm[i], c) for i, c in w]
                    key = _reanchored_key(g) if g else ()
                    if best is None or key < best:
                        best = key
                canon.add(best)
    return len(canon)


# ---------------------------------------------------------------------------
# The 8-site box B_1 = {0,1}^3 (half-integer centre): same order-48 group,
# small enough for full brute force; machine-check data for the test.
# ---------------------------------------------------------------------------

B1_SITES: tuple[tuple[int, int, int], ...] = tuple(itertools.product(range(2), repeat=3))
B1_INDEX = {s: i for i, s in enumerate(B1_SITES)}
B1_CTR = tuple((s[0] * 2 - 1, s[1] * 2 - 1, s[2] * 2 - 1) for s in B1_SITES)


def _b1_perm_table(p: tuple[int, ...], signs: tuple[int, ...]) -> tuple[int, ...]:
    table = []
    for idx in range(8):
        c = B1_CTR[idx]
        img = tuple(signs[i] * c[p[i]] for i in range(3))
        table.append(B1_INDEX[tuple((v + 1) // 2 for v in img)])
    return tuple(table)


B1_PERMS = tuple(_b1_perm_table(p, s) for p, s in GROUP_SIGNS)


def b1_column_orbit_counts(smax: int) -> list[int]:
    fixes = []
    for perm in B1_PERMS:
        seen = [False] * 8
        cycles = []
        for i in range(8):
            if seen[i]:
                continue
            j, length = i, 0
            while not seen[j]:
                seen[j] = True
                j = perm[j]
                length += 1
            cycles.append(length)
        fixes.append(fix_counts_from_cycles(cycles, smax))
    orbits = []
    for s in range(smax + 1):
        total = sum(f[s] for f in fixes)
        if total % GROUP_ORDER:
            raise AssertionError(f"B1 column Burnside not integral at s={s}")
        orbits.append(total // GROUP_ORDER)
    return orbits


def b1_brute_column_orbits(smax: int) -> int:
    canon: set = set()
    for k in range(smax + 1):
        for sites in itertools.combinations(range(8), k):
            for codes in itertools.product((1, 2, 3), repeat=k):
                w = tuple(zip(sites, codes))
                best = None
                for perm in B1_PERMS:
                    key = tuple(sorted((perm[i], c) for i, c in w))
                    if best is None or key < best:
                        best = key
                canon.add(best)
    return len(canon)


def b1_row_orbit_counts(smax: int) -> list[int]:
    """Affine DP on doubled centred coordinates: A(D) = g(D) - 2t, t in
    {-1,0,1}^3; anchored means touching D_i = -1 on every axis."""
    per_g: list[list[int]] = []
    for p, signs in GROUP_SIGNS:
        fix = [1] * (smax + 1)
        for t in itertools.product((-1, 0, 1), repeat=3):
            t2 = (2 * t[0], 2 * t[1], 2 * t[2])
            cycles, seen_cycles, done = [], set(), set()
            for start in range(8):
                if start in done:
                    continue
                walk, posmap, cur = [start], {start: 0}, start
                while True:
                    cy = B1_CTR[cur]
                    ny = tuple(signs[i] * cy[p[i]] - t2[i] for i in range(3))
                    if not all(v in (-1, 1) for v in ny):
                        break
                    nxt = B1_INDEX[tuple((v + 1) // 2 for v in ny)]
                    if nxt in posmap:
                        cyc = frozenset(walk[posmap[nxt]:])
                        if cyc not in seen_cycles:
                            seen_cycles.add(cyc)
                            cycles.append(cyc)
                        break
                    posmap[nxt] = len(walk)
                    walk.append(nxt)
                    cur = nxt
                done.update(walk)
            dp = [[0] * 8 for _ in range(smax + 1)]
            dp[0][0] = 1
            for cyc in cycles:
                length = len(cyc)
                cov = 0
                for i in cyc:
                    c = B1_CTR[i]
                    if c[0] == -1:
                        cov |= 1
                    if c[1] == -1:
                        cov |= 2
                    if c[2] == -1:
                        cov |= 4
                new = [row[:] for row in dp]
                for n in range(length, smax + 1):
                    base = dp[n - length]
                    for c0 in range(8):
                        if base[c0]:
                            new[n][c0 | cov] += 3 * base[c0]
                dp = new
            for s in range(smax + 1):
                fix[s] += sum(dp[n][7] for n in range(1, s + 1))
        per_g.append(fix)
    orbits = []
    for s in range(smax + 1):
        total = sum(f[s] for f in per_g)
        if total % GROUP_ORDER:
            raise AssertionError(f"B1 row Burnside not integral at s={s}")
        orbits.append(total // GROUP_ORDER)
    return orbits


def b1_brute_row_orbits(smax: int) -> int:
    def rekey(pairs):
        pts = [(B1_SITES[i], c) for i, c in pairs]
        shift = [min(p[0][a] for p in pts) for a in range(3)]
        return tuple(sorted(
            ((p[0][0] - shift[0], p[0][1] - shift[1], p[0][2] - shift[2]), p[1])
            for p in pts))

    canon = set()
    for k in range(0, smax + 1):
        for sites in itertools.combinations(range(8), k):
            if k:
                pts = [B1_SITES[i] for i in sites]
                if not (min(p[0] for p in pts) == 0
                        and min(p[1] for p in pts) == 0
                        and min(p[2] for p in pts) == 0):
                    continue
            for codes in itertools.product((1, 2, 3), repeat=k):
                w = tuple(zip(sites, codes))
                best = None
                for perm in B1_PERMS:
                    g = [(perm[i], c) for i, c in w]
                    key = rekey(g) if g else ()
                    if best is None or key < best:
                        best = key
                canon.add(best)
    return len(canon)


# ---------------------------------------------------------------------------
# Orbit-representative column stream: canonical word = orbit minimum under
# the lexicographic (support mask, letters) order on the 27 base sites.
# ---------------------------------------------------------------------------

def _apply_perm_mask(mask: int, perm: tuple[int, ...]) -> int:
    out = 0
    while mask:
        bit = mask & -mask
        mask ^= bit
        out |= 1 << perm[bit.bit_length() - 1]
    return out


def canonical_supports(k: int) -> Iterator[tuple[tuple[int, ...], list[tuple[int, ...]]]]:
    """Yield (sorted sites, stabiliser site-permutation tables) for every
    canonical (orbit-minimum mask) k-subset of the 27 sites."""
    seen: set[int] = set()
    for combo in itertools.combinations(range(27), k):
        mask = 0
        for i in combo:
            mask |= 1 << i
        if mask in seen:
            continue
        orbit = [_apply_perm_mask(mask, p) for p in PERMS27]
        seen.update(orbit)
        canonical = min(orbit)
        sites = []
        j = canonical
        pos = 0
        while j:
            if j & 1:
                sites.append(pos)
            j >>= 1
            pos += 1
        stab = [p for p in PERMS27 if _apply_perm_mask(canonical, p) == canonical]
        yield tuple(sites), stab


def iter_representative_columns(smax: int) -> Iterator[tuple[int, int]]:
    """Stream one canonical representative per G-orbit of class words.

    The canonical word is the orbit minimum under the order that compares
    the 27-bit support mask first and the base-4 letter packing (letters of
    the support sites in increasing site order) second.  A word is
    canonical iff its support mask is the orbit-minimum mask and its
    letters are minimal over the support stabiliser.  The identity word is
    its own orbit and is yielded first.
    """
    yield (0, 0)
    for size in range(1, smax + 1):
        for sites, stab in canonical_supports(size):
            k = len(sites)
            bits = [BASE_BITS[i] for i in sites]
            if len(stab) == 1:
                for codes in itertools.product((1, 2, 3), repeat=k):
                    x = z = 0
                    for b, code in zip(bits, codes):
                        if code & 1:
                            x |= b
                        if code & 2:
                            z |= b
                    yield (x, z)
                continue
            shifts = [k - 1 - j for j in range(k)]
            # For each stabiliser element: position permutation of letters.
            posmaps = []
            site_pos = {s: j for j, s in enumerate(sites)}
            for p in stab:
                posmaps.append(tuple(site_pos[p[s]] for s in sites))
            for codes in itertools.product((1, 2, 3), repeat=k):
                key0 = 0
                for j in range(k):
                    key0 |= codes[j] << (2 * shifts[j])
                ok = True
                for pm in posmaps:
                    keyg = 0
                    for j in range(k):
                        keyg |= codes[j] << (2 * shifts[pm[j]])
                    if keyg < key0:
                        ok = False
                        break
                if ok:
                    x = z = 0
                    for b, code in zip(bits, codes):
                        if code & 1:
                            x |= b
                        if code & 2:
                            z |= b
                    yield (x, z)


def representative_column_count(smax: int) -> int:
    return sum(1 for _ in iter_representative_columns(smax))


# ---------------------------------------------------------------------------
# Column echelon over F_p (wave-14 engine, generalised to any column stream).
# ---------------------------------------------------------------------------

class LeadTable:
    """Fixed-capacity linear-probing map key -> pivot index."""

    __slots__ = ("cap", "mask", "keys", "vals")

    def __init__(self, cap_bits: int) -> None:
        self.cap = 1 << cap_bits
        self.mask = self.cap - 1
        self.keys = array("q", bytes(8 * self.cap))
        self.vals = array("i", bytes(4 * self.cap))

    def slot_of(self, key: int) -> int:
        slot = (key * 0x9E3779B97F4A7C15) & 0x1FFFFFFFFFFFFF & self.mask
        keys = self.keys
        while True:
            k = keys[slot]
            if k == key or k == 0:
                return slot
            slot = (slot + 1) & self.mask




def eliminate(stage: str, smax: int, columns: Iterator[tuple[int, int]],
              total_columns: int, prime: int, budget_s: float,
              cap_bits: int, guard_pivots: int,
              log_interval: int = 100_000) -> dict[str, object]:
    """Stream the given columns; return the exact F_p rank and measurements.

    Identical certificate engine to wave 14 (static descending
    (support size, packed blob) lead priority; pivot count = F_p column
    rank, independent of lead order).  Raises ResourceWall on an OBSERVED
    process_time or RSS breach, carrying exact numbers.
    """
    t0 = process_time()
    gc.disable()
    lead = LeadTable(cap_bits)
    offsets = array("q", [0])
    keys = array("q")
    coeffs = array("i")
    digest = hashlib.sha256()
    pcount = stored = maxlen = reductions = ops = 0
    cid = -1
    try:
        for x, z in columns:
            cid += 1
            if not x and not z:
                continue  # identity column: zero image
            w: dict[int, int] = {}
            for (nx, nz), c in image_terms(x, z).items():
                assert not (nx | nz) & ~EXT_MASK, "Lemma A violated: word left C_R"
                k = anchor_key(nx, nz)
                assert (k >> 60) <= smax + 1, "Lemma A violated: support size > s+1"
                w[k] = (w.get(k, 0) + c) % prime
            w = {k: v for k, v in w.items() if v}
            while w:
                lk = max(w)
                slot = lead.slot_of(lk)
                if lead.keys[slot] == 0:
                    inv = pow(w[lk], prime - 2, prime)
                    del w[lk]
                    for k, v in w.items():
                        keys.append(k)
                        coeffs.append(v * inv % prime)
                    offsets.append(len(keys))
                    lead.keys[slot] = lk
                    lead.vals[slot] = pcount
                    digest.update(cid.to_bytes(8, "little"))
                    digest.update(lk.to_bytes(16, "little"))
                    pcount += 1
                    stored += len(w)
                    if len(w) > maxlen:
                        maxlen = len(w)
                    if pcount > guard_pivots:
                        raise ResourceWall(
                            stage,
                            f"pivot count {pcount} exceeded plan guard {guard_pivots}",
                            cid + 1, total_columns, process_time() - t0, peak_rss_bytes())
                    break
                pv = lead.vals[slot]
                reductions += 1
                sc = w.pop(lk)
                lo = offsets[pv]
                hi = offsets[pv + 1]
                ops += hi - lo
                g = w.get
                for i in range(lo, hi):
                    k = keys[i]
                    nv = (g(k, 0) - sc * coeffs[i]) % prime
                    if nv:
                        w[k] = nv
                    else:
                        w.pop(k, None)
            if cid % log_interval == 0:
                el = process_time() - t0
                rss = peak_rss_bytes()
                print(f"  [{stage} p{prime % 1000}] col {cid:>10,} piv={pcount:>10,} "
                      f"nnz={stored:>12,} avg={stored / max(pcount, 1):5.2f} "
                      f"red={reductions:>10,} t={el:7.0f}s rss={rss / 2**20:6.0f}MiB",
                      flush=True)
                if el > budget_s:
                    raise ResourceWall(
                        stage, f"process_time budget {budget_s}s exceeded",
                        cid + 1, total_columns, el, rss)
                if rss > RSS_WALL_BYTES:
                    raise ResourceWall(
                        stage, f"peak RSS {rss / 2**20:.1f} MiB exceeded the 8 GiB wall",
                        cid + 1, total_columns, el, rss)
    finally:
        gc.enable()
    el = process_time() - t0
    rss = peak_rss_bytes()
    if el > budget_s:
        raise ResourceWall(stage, f"process_time budget {budget_s}s exceeded at stage end",
                           cid + 1, total_columns, el, rss)
    return {
        "claim_tag": "[COMPUTATION]",
        "stage": stage,
        "smax": smax,
        "prime": prime,
        "columns": cid + 1,
        "columns_expected": total_columns,
        "rank_Fp": pcount,
        "pivot_vectors_nnz_excluding_leads": stored,
        "average_pivot_length": round(stored / max(pcount, 1), 6),
        "maximum_pivot_length": maxlen,
        "reductions": reductions,
        "csr_ops": ops,
        "elapsed_seconds_process_time": round(el, 6),
        "budget_seconds_process_time": budget_s,
        "peak_rss_bytes": rss,
        "peak_rss_mib": round(rss / 1024**2, 3),
        "pivot_trace_sha256_of_column_lead_stream": digest.hexdigest(),
        "lead_priority": ("static descending (support size, packed blob); the "
                          "pivot count of a column echelon is lead-order independent"),
    }


# ---------------------------------------------------------------------------
# Defect traps (cheap, decisive, run before anything expensive).
# ---------------------------------------------------------------------------

def defect_traps() -> dict[str, object]:
    origin = 1 << EXT_INDEX[(0, 0, 0)]
    dx = image_terms(origin, 0)
    dz = image_terms(0, origin)
    # h = X_0 + sum_i Z_0 Z_{e_i}: pi D h must vanish identically.
    summed: dict[int, int] = {}
    h_terms = [(origin, 0)]
    for axis in range(3):
        nb = [0, 0, 0]
        nb[axis] = 1
        h_terms.append((0, origin | (1 << EXT_INDEX[tuple(nb)])))
    for hx, hz in h_terms:
        for (nx, nz), c in image_terms(hx, hz).items():
            k = anchor_key(nx, nz)
            summed[k] = summed.get(k, 0) + c
    summed = {k: v for k, v in summed.items() if v}
    return {
        "claim_tag": "[COMPUTATION]",
        "D_of_X0_terms": sorted(dx.values()),
        "D_of_X0_count": len(dx),
        "D_of_Z0": {f"{k[0]:b}_{k[1]:b}": v for k, v in dz.items()},
        "pi_D_h_empty": not summed,
        "pi_D_h_residual_terms": len(summed),
        "h_word_support_sizes": sorted(bin(x | z).count("1") for x, z in h_terms),
    }


# ---------------------------------------------------------------------------
# Orchestration.
# ---------------------------------------------------------------------------

def check(name: str, passed: bool, detail: str) -> dict[str, object]:
    return {"name": name, "passed": bool(passed), "detail": detail}


def main() -> int:
    started = process_time()
    checks: list[dict[str, object]] = []
    print("e144_octahedral: order-48 point-group feasibility + execution", flush=True)

    # ----- FEASIBILITY (exact integer arithmetic, before any big run) -------
    print("[preflight] Burnside arithmetic ...", flush=True)
    t0 = process_time()
    columns_by_size = [math.comb(27, k) * 3**k for k in range(7)]
    cols6 = sum(columns_by_size)
    cols5 = sum(columns_by_size[:6])
    rank_S = [rank_S_closed_form(s) for s in range(7)]
    col_orbits, col_fixes = column_orbit_counts(6)
    row_orbits, row_fixes = row_orbit_counts(6)
    if process_time() - t0 > PREFLIGHT_BUDGET_S:
        print("WARNING: preflight arithmetic exceeded budget", flush=True)

    checks.append(check(
        "group_order_48_distinct_site_actions",
        len(set(PERMS27)) == GROUP_ORDER,
        f"{len(set(PERMS27))} distinct site permutations on the 27 sites"))
    checks.append(check(
        "column_burnside_integral_all_s",
        all(sum(f[s] for f in col_fixes) % GROUP_ORDER == 0 for s in range(7)),
        "sum of fixes divisible by 48 at every s<=6"))
    checks.append(check(
        "row_burnside_integral_all_s",
        all(sum(f[s] for f in row_fixes) % GROUP_ORDER == 0 for s in range(7)),
        "sum of fixes divisible by 48 at every s<=6"))
    # (e, t=0) identity check: the affine DP must reproduce the rank-S rows.
    e_elem = next((p, s) for p, s in GROUP_SIGNS
                  if p == (0, 1, 2) and s == (1, 1, 1))
    ok_identity = True
    for s in range(7):
        dp = row_fix_from_cycles(affine_in_box_cycles(*e_elem, (0, 0, 0)), s)
        if dp[s] != rank_S[s] - 1:
            ok_identity = False
    checks.append(check(
        "identity_affine_dp_equals_rankS_minus_identity_all_s", ok_identity,
        "at (g=e, t=0) the row DP returns 1+sum_k A_k 3^k rows for every s<=6"))

    print("[preflight] brute-force validations ...", flush=True)
    bf_cols3 = brute_column_orbits(3)
    bf_rows3 = brute_row_orbits(3)
    checks.append(check(
        "brute_column_orbits_agree_s3", bf_cols3 == col_orbits[3],
        f"brute {bf_cols3:,} vs Burnside {col_orbits[3]:,}"))
    checks.append(check(
        "brute_row_orbits_agree_s3", bf_rows3 == row_orbits[3],
        f"brute {bf_rows3:,} vs Burnside {row_orbits[3]:,}"))
    b1_cols = b1_column_orbit_counts(8)
    b1_cols_bf = b1_brute_column_orbits(8)
    b1_rows = b1_row_orbit_counts(4)
    b1_rows_bf = b1_brute_row_orbits(4)
    checks.append(check(
        "b1_column_burnside_vs_bruteforce_s8", b1_cols[8] == b1_cols_bf,
        f"B1 brute {b1_cols_bf:,} vs Burnside {b1_cols[8]:,}"))
    checks.append(check(
        "b1_row_burnside_vs_bruteforce_s4", b1_rows[4] == b1_rows_bf,
        f"B1 brute {b1_rows_bf:,} vs Burnside {b1_rows[4]:,}"))

    checks.append(check(
        "closed_forms_match_frozen_anchors",
        class_column_count(4) == CERTIFIED_S4_COLUMNS
        and rank_S_closed_form(4) == CERTIFIED_S4_RANK_S
        and class_column_count(5) == CERTIFIED_S5_COLUMNS
        and rank_S_closed_form(5) == CERTIFIED_S5_RANK_S,
        "s<=4/s<=5 column counts and rank S reproduce the frozen wave-13/14 values"))
    traps = defect_traps()
    checks.append(check(
        "defect_traps",
        traps["D_of_X0_terms"] == [-1] * 6
        and traps["D_of_X0_count"] == 6
        and traps["pi_D_h_empty"]
        and traps["h_word_support_sizes"] == [1, 2, 2, 2],
        "D(X_0) has six bond terms all -1; pi D h = 0 exactly on orbit sums"))

    # ----- EXECUTION: full-class regressions -------------------------------
    stages: dict[str, dict[str, object]] = {}
    walls: list[dict[str, object]] = []

    def run_stage(stage_key: str, smax: int, kind: str) -> None:
        cap_bits, guard, budget, log_interval = STAGE_PLAN[stage_key]
        if kind == "full":
            total = class_column_count(smax)
        else:
            total = col_orbits[smax]
        runs = []
        for prime in PRIMES:
            print(f"[{stage_key}] smax={smax} kind={kind} cols={total:,} "
                  f"prime={prime}", flush=True)
            gen = (iter_class_columns(smax) if kind == "full"
                   else iter_representative_columns(smax))
            try:
                runs.append(eliminate(stage_key, smax, gen, total, prime,
                                      budget, cap_bits, guard, log_interval))
            except ResourceWall as wall:
                walls.append(wall.record)
                runs.append({"claim_tag": "[COMPUTATION]", "walled": True,
                             "stage": stage_key, "prime": prime,
                             **wall.record})
        stages[stage_key] = {
            "claim_tag": "[COMPUTATION]",
            "kind": kind,
            "smax": smax,
            "columns": total,
            "runs": runs,
        }

    run_stage("full_s2", 2, "full")
    run_stage("full_s3", 3, "full")
    run_stage("full_s4", 4, "full")

    s4_runs = stages["full_s4"]["runs"]
    s4_ranks = [r.get("rank_Fp") for r in s4_runs if not r.get("walled")]
    checks.append(check(
        "s4_regression_reproduces_frozen_rank",
        len(s4_ranks) == 2 and s4_ranks[0] == CERTIFIED_S4_RANK_M
        and s4_ranks[1] == CERTIFIED_S4_RANK_M,
        f"full s<=4 echelon ranks {s4_ranks} vs frozen {CERTIFIED_S4_RANK_M}"))
    s3_ranks = [r.get("rank_Fp") for r in stages["full_s3"]["runs"]
                if not r.get("walled")]
    s2_ranks = [r.get("rank_Fp") for r in stages["full_s2"]["runs"]
                if not r.get("walled")]

    # ----- EXECUTION: orbit-representative echelons -------------------------
    print(f"[rep_s5] columns={col_orbits[5]:,} (full-class target "
          f"{CERTIFIED_S5_RANK_M:,})", flush=True)
    run_stage("rep_s5", 5, "rep")
    print(f"[rep_s6] columns={col_orbits[6]:,} (target rank {rank_S[6] - 2:,})",
          flush=True)
    run_stage("rep_s6", 6, "rep")

    rep_counts_ok5 = stages["rep_s5"]["runs"] and all(
        r.get("columns") == col_orbits[5] for r in stages["rep_s5"]["runs"]
        if not r.get("walled"))
    rep_counts_ok6 = stages["rep_s6"]["runs"] and all(
        r.get("columns") == col_orbits[6] for r in stages["rep_s6"]["runs"]
        if not r.get("walled"))
    checks.append(check(
        "representative_streams_match_burnside_counts",
        bool(rep_counts_ok5 and rep_counts_ok6),
        f"rep streams emitted {col_orbits[5]:,}/{col_orbits[6]:,} columns = "
        "the Burnside orbit counts (canonical-form enumeration is exact)"))

    rep5_ranks = [r.get("rank_Fp") for r in stages["rep_s5"]["runs"]
                  if not r.get("walled")]
    rep6_ranks = [r.get("rank_Fp") for r in stages["rep_s6"]["runs"]
                  if not r.get("walled")]
    checks.append(check(
        "rep_ranks_within_orbit_caps",
        all(r <= col_orbits[5] for r in rep5_ranks)
        and all(r <= col_orbits[6] for r in rep6_ranks),
        f"rep s<=5 ranks {rep5_ranks} <= {col_orbits[5]:,}; "
        f"rep s<=6 ranks {rep6_ranks} <= {col_orbits[6]:,}"))

    # ----- DECISION: cap lemma + memory/time models -------------------------
    required_pivots = rank_S[6] - 2
    cap_holds = col_orbits[6] < required_pivots
    checks.append(check(
        "cap_lemma_blocks_representative_decision", cap_holds,
        f"{col_orbits[6]:,} orbit representatives < {required_pivots:,} "
        "required pivots: a representative echelon certifies at most "
        f"{100 * col_orbits[6] // required_pivots}% of the target rank"))

    lead_cap_bits_6 = 28  # 268,435,456 slots; load 191,949,608/2^28 = 0.715
    lead_bytes = 12 << lead_cap_bits_6
    csr_bytes = required_pivots * FROZEN_S5_PIVOT_NNZ * 12 // FROZEN_S5_PIVOTS
    whole_process_bytes = (required_pivots * FROZEN_S5_PEAK_RSS
                           // FROZEN_S5_PIVOTS)
    total_model_bytes = lead_bytes + csr_bytes
    fits = total_model_bytes <= RSS_WALL_BYTES and whole_process_bytes <= RSS_WALL_BYTES
    checks.append(check(
        "full_s6_memory_models_exceed_8gib", not fits,
        f"lead+CSR model {total_model_bytes:,} B and whole-process model "
        f"{whole_process_bytes:,} B both exceed the 8 GiB = "
        f"{RSS_WALL_BYTES:,} B gate"))
    rate = CERTIFIED_S5_COLUMNS / FROZEN_S5_SECONDS  # columns per second
    projected_seconds_per_prime = int(cols6 / rate)

    # Certified consequences of the representative echelons.
    bound6 = min(rep6_ranks) if rep6_ranks else None
    quotient_upper = (rank_S[6] - bound6 - 1) if bound6 is not None else None
    decided = (bound6 is not None and bound6 == required_pivots)
    checks.append(check(
        "s6_quotient_not_decided_by_orbit_reduction", not decided,
        "the s<=6 quotient remains bracketed: 1 <= Q(C(3,2,6)) <= "
        f"{quotient_upper:,}" if quotient_upper is not None else
        "no representative run completed"))

    verdict = ("decided_quotient_1" if decided
               else "wall_certificate_s6_undecided_by_octahedral_reduction")

    preflight_seconds = process_time() - started

    payload_data = {
        "claim_tags": ["[COMPUTATION]", "[THEOREM-elementary]"],
        "setting": {
            "class": "C(3,2,s): rational span of ordered-Pauli words supported "
                     "in B_2 = {0,1,2}^3 with word support size <= s",
            "group": "octahedral point group O_h (48 signed coordinate "
                     "permutations about the box centre), the Z^3 lattice "
                     "point group preserving the 3x3x3 box",
            "hamiltonian": "H = sum X_x + sum_<xy> Z_x Z_y on Z^3; "
                           "D = (1/2)[H, .]; pi = translation-orbit sum",
            "quotient_identity": "Q(C) = rank S_C - rank M_C - 1 "
                                 "(wave-12/13 framework)",
        },
        "feasibility": {
            "claim_tag": "[COMPUTATION]",
            "column_counts_by_size": columns_by_size,
            "columns_s6": cols6,
            "columns_s5": cols5,
            "rank_S_by_s": rank_S,
            "column_burnside": {
                "claim_tag": "[THEOREM-elementary]+[COMPUTATION]",
                "group_order": GROUP_ORDER,
                "orbit_counts_by_s": col_orbits,
                "per_element_fix_sum_by_s": [sum(f[s] for f in col_fixes)
                                             for s in range(7)],
                "fix_by_element": col_fixes,
                "reduction_factor_s6": f"{cols6}/{col_orbits[6]}",
                "brute_force_validation_s3": bf_cols3,
                "brute_force_validation_s3_agrees": bf_cols3 == col_orbits[3],
            },
            "row_burnside": {
                "claim_tag": "[THEOREM-elementary]+[COMPUTATION]",
                "orbit_counts_by_s": row_orbits,
                "per_element_fix_sum_by_s": [sum(f[s] for f in row_fixes)
                                             for s in range(7)],
                "fix_by_element": row_fixes,
                "rows_by_s": rank_S,
                "reduction_factor_s6": f"{rank_S[6]}/{row_orbits[6]}",
                "identity_element_check": "at (g=e, t=0) the affine DP "
                                          "reproduces 1+sum_k A_k 3^k for every s",
                "brute_force_validation_s3": bf_rows3,
                "brute_force_validation_s3_agrees": bf_rows3 == row_orbits[3],
            },
            "b1_box_machine_check": {
                "claim_tag": "[COMPUTATION]",
                "box": "B_1 = {0,1}^3 (8 sites, half-integer centre)",
                "column_orbits_by_s": b1_cols,
                "column_orbits_bruteforce_s8": b1_cols_bf,
                "row_orbits_by_s": b1_rows,
                "row_orbits_bruteforce_s4": b1_rows_bf,
            },
            "defect_traps": traps,
        },
        "regression_full_class": {
            "s2": stages["full_s2"], "s3": stages["full_s3"],
            "s4_frozen_anchor": stages["full_s4"],
            "s4_frozen_rank_M": CERTIFIED_S4_RANK_M,
        },
        "rep_calibration_s5": {
            **stages["rep_s5"],
            "full_class_target_rank": CERTIFIED_S5_RANK_M,
            "full_class_columns": CERTIFIED_S5_COLUMNS,
            "stall_fraction_of_target": (
                f"{min(rep5_ranks) if rep5_ranks else 0}/{CERTIFIED_S5_RANK_M}"
                if rep5_ranks else "walled"),
        },
        "headline_s6_reps": {
            **stages["rep_s6"],
            "required_pivots_for_decision": required_pivots,
            "certified_lower_bound_rank_Q_M6": bound6,
            "certified_quotient_bracket": (
                [1, quotient_upper] if quotient_upper is not None else None),
            "bracket_statement": (
                "1 <= Q(C(3,2,6)) <= rank S_6 - min_p rank_Fp(M|reps) - 1; the "
                "lower end is the Hamiltonian density h (analytic, wave-14 "
                "argument), the upper end follows from rank_Q M >= rank_Fp of "
                "any integer submatrix"),
        },
        "decision": {
            "claim_tag": "[COMPUTATION]",
            "required_pivots": required_pivots,
            "cap_lemma": {
                "claim_tag": "[THEOREM-elementary]",
                "statement": "a column echelon fed only G-orbit "
                             "representatives certifies at most #orbits pivots; "
                             "#orbits < rank S_6 - 2, so orbit-representative "
                             "compression cannot decide the s<=6 quotient",
                "orbit_representatives": col_orbits[6],
                "fraction_of_target": f"{col_orbits[6]}/{required_pivots}",
                "fraction_percent_floor": 100 * col_orbits[6] // required_pivots,
            },
            "memory_models_bytes": {
                "claim_tag": "[COMPUTATION]",
                "lead_table": lead_bytes,
                "lead_table_cap_bits": lead_cap_bits_6,
                "csr_pivot_store": csr_bytes,
                "csr_model_note": "required pivots x frozen s<=5 measured "
                                  "average pivot nnz (21.248964 excluding lead) "
                                  "x 12 bytes/key+coefficient",
                "total_lead_plus_csr": total_model_bytes,
                "whole_process_model": whole_process_bytes,
                "whole_process_note": "required pivots x frozen s<=5 peak RSS "
                                      "(2,404,089,856 B) / frozen pivots",
                "wall_bytes": RSS_WALL_BYTES,
                "fits_8gib": fits,
            },
            "time_model": {
                "claim_tag": "[COMPUTATION]",
                "frozen_s5_columns_per_second": rate,
                "projected_seconds_per_prime": projected_seconds_per_prime,
                "projected_hours_both_primes": round(
                    2 * projected_seconds_per_prime / 3600, 2),
            },
            "surviving_fraction_analysis": {
                "claim_tag": "[COMPUTATION]+[CONJECTURE]",
                "columns_surviving_as_representatives": f"{col_orbits[6]}/{cols6}",
                "rank_target_reduction": "none: the image of M_6 is G-stable of "
                                         "dimension rank S_6 - 2, which no "
                                         "sector/orbit bookkeeping shrinks; the "
                                         "certificate still needs "
                                         f"{required_pivots:,} pivots",
                "rows_surviving_as_orbits": f"{row_orbits[6]}/{rank_S[6]}",
                "conjectural_equivariant_route": (
                    "a G-equivariant elimination storing one pivot row per "
                    "ROW orbit (~"
                    f"{row_orbits[6] * FROZEN_S5_PIVOT_NNZ * 12 // FROZEN_S5_PIVOTS:,}"
                    " B CSR) plus a lead structure over lead orbits would "
                    "plausibly fit the 8 GiB gate, but requires on-the-fly "
                    "exact reconstruction of stored pivot rows under every "
                    "group element with a G-compatible pivot order; NOT "
                    "implemented here, recorded as the identified follow-up"),
            },
            "verdict": verdict,
        },
        "unresolved": [
            "the s<=6 class quotient Q(C(3,2,6)) is bracketed in "
            + (f"[1, {quotient_upper}]" if quotient_upper is not None else "[1, ?]")
            + " and NOT decided; the full radius-2 Z^3 box remains untouched",
            "no all-size charge theorem follows from these finite statements",
        ],
        "walls_observed": walls,
    }

    provenance = {
        "generated_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "script": SCRIPT,
        "interpreter": ".venv/bin/python",
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "method": ("exact integer Burnside arithmetic (cycle DP for columns; "
                   "affine A_t cycle DP for anchored rows); canonical-form "
                   "orbit representatives; wave-14 column echelon over F_p at "
                   "two primes; decision by cap lemma + frozen-calibrated "
                   "memory models"),
        "certifying_arithmetic": "Python integers; finite-field elimination "
                                 "modulo the listed primes; no floats in any "
                                 "certified quantity",
        "clock": "time.process_time() only; never wall clock or signal alarms",
        "rss_measurement": "process-lifetime ru_maxrss at measurement end; "
                           "conservative",
        "budgets_seconds_process_time": {
            "preflight": PREFLIGHT_BUDGET_S,
            "full_s2_per_prime": S2_PRIME_BUDGET_S,
            "full_s3_per_prime": S3_PRIME_BUDGET_S,
            "full_s4_per_prime": S4_PRIME_BUDGET_S,
            "rep_s5_per_prime": REP5_PRIME_BUDGET_S,
            "rep_s6_per_prime": REP6_PRIME_BUDGET_S,
        },
        "rss_wall_bytes": RSS_WALL_BYTES,
        "frozen_anchors": {
            "s4": "results/integrability/extensive_3d_r2.json",
            "s5": "results/integrability/extensive_3d_s5.json",
        },
        "preflight_seconds_process_time": round(preflight_seconds, 3),
        "total_elapsed_seconds_process_time": round(process_time() - started, 3),
    }

    RESULT.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(
        {"checks": checks, "data": payload_data, "provenance": provenance},
        indent=2) + "\n")

    failed = [c for c in checks if not c["passed"]]
    for c in checks:
        status = "PASS" if c["passed"] else "FAIL"
        print(f"[{status}] {c['name']}: {c['detail']}", flush=True)
    print(f"verdict: {verdict}", flush=True)
    print(f"total process_time {process_time() - started:.1f}s; "
          f"{len(failed)} failing checks", flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
