#!/usr/bin/env python3
"""Column-oriented exact elimination for the Z^3 R=2 support-size class s<=5.

Wave 13 (`experiments/e125_extensive_3d_r2.py`,
`proofs/extensive_3d_r2.md`) certified the class quotient

    dim ker(pi o D|_C) / (C cap ker pi + Q I) = rank S_C - rank M_C - 1

for the support-size classes C(3,2,s), s<=4, and hit an OBSERVED 5 GiB RSS
wall probing s<=5: the row-oriented streamed-CSR enumeration stored every
column image and died at 1,654,784 of 21,121,156 columns after 171.8 s at
5,133 MiB process-lifetime peak RSS.  That baseline is what this script
must beat.

Method: *column-oriented* exact elimination over F_p, calibrated on the
certified s<=4 class (it reproduces rank 822,332 exactly):

* columns (class words) are STREAMED one at a time and never stored; only
  pivot (echelon) vectors are kept, in a compact CSR store (array('q')
  keys, array('i') coefficients, array('q') offsets);
* each pivot vector EXCLUDES its own lead entry (recoverable from the lead
  map), so a reduction updates the working vector in place with no copy;
* the lead map is an open-addressing hash table over two flat arrays
  (multiplicative hashing, fixed capacity, no per-entry Python objects);
* row keys are single integers under a PROVABLY INJECTIVE packing:
  digit = ((dx*5+dy)*5+dz)*4 + (code-1) < 500 with dx,dy,dz in 0..4 and
  code in 1..3, at most smax+1 = 6 sorted digits, so the packed part is
  < 500^6 < 2^54 and key = size<<54 | packed is injective.  (The wave-14
  calibration scratch used a 44-bit size field, which is NOT injective for
  sizes >= 5; this script's (size, packed) ordering is identical but the
  packing is collision-free.);
* lead priority is the static descending (support size, packed blob) order
  of these keys (max of the working dict).  The pivot count of a column
  echelon over a field is the rank and is independent of the lead order;
  the priority only controls fill-in (predecessor calibration: average
  pivot length 18.2 including its normalized lead on s<=4).

Certificate structure (identical to wave 13, no new mathematics): with
rank S_C computed exactly by the anchored-orbit closed form and the
analytic kernel {I, h} giving rank_Q M_C <= rank S_C - 2, a modular rank
equal to rank S_C - 2 at two primes forces rank_Q M_C = rank S_C - 2 by
the sandwich rank_Fp <= rank_Q, certifying quotient = 1: modulo trivial
translation shifts and the identity, the only class density with vanishing
projected commutator is the Hamiltonian density h.

Mandatory preflight: the s<=4 certified values (columns 1,503,766,
rank S 822,334, rank M 822,332 at both primes) are re-verified as a
regression BEFORE the new class runs.  A static greedy row-discovery
lower bound (certified correct by a triangular-minor argument) is also
recorded: it stalls strictly below the exact rank on the calibration
class, so the exact elimination is the certificate, not the greedy count.

Budgets (predeclared, process_time only, never wall clock; RSS under
6 GB = 6e9 bytes):
  * closed-form arithmetic + defect traps:                    600 s
  * s<=4 regression, 2,400 s total (1,200 s per prime);
  * s<=4 greedy lower bounds, 900 s total (450 s per rule);
  * s<=5 headline, 7,200 s per prime.
Any budget or RSS breach raises ResourceWall and is recorded as an
OBSERVED wall with exact numbers; no rank is claimed for a walled case.

Run: PYTHONPATH=src .venv/bin/python experiments/e135_extensive_3d_s5.py
"""
from __future__ import annotations

import gc
import hashlib
import itertools
import json
import math
import platform
import resource
import time
from array import array
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterator

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "integrability" / "extensive_3d_s5.json"
PRIOR_RESULT = ROOT / "results" / "integrability" / "extensive_3d_r2.json"
SCRIPT = "experiments/e135_extensive_3d_s5.py"

PRIMES = (2_147_483_647, 2_147_483_629)

# Predeclared budgets: process_time seconds (never wall clock), RSS bytes.
RSS_WALL_BYTES = 6_000_000_000  # 6 GB decimal, stricter than 6 GiB
PREFLIGHT_BUDGET_S = 600.0
REGRESSION_BUDGET_S = 2_400.0  # total, both primes
GREEDY_BUDGET_S = 900.0  # total, both rules
S5_PRIME_BUDGET_S = 7_200.0  # per prime
RSS_CHECK_COLUMN_INTERVAL = 1_000_000

# Frozen wave-13 certified anchors for the regression stage.
CERTIFIED_S4_COLUMNS = 1_503_766
CERTIFIED_S4_RANK_S = 822_334
CERTIFIED_S4_RANK_M = 822_332

# Closed forms for the new class (re-derived and cross-checked at runtime).
S5_COLUMNS = 21_121_156
S5_RANK_S = 14_757_412
S5_TARGET_RANK = S5_RANK_S - 2  # analytic {I,h} kernel upper bound = 14,757,410

# Lead-table capacities (slots); guards leave substantial unused capacity.
CAP_BITS_S4 = 21  # 2,097,152 slots vs 1,470,000 guard
CAP_BITS_S5 = 25  # 33,554,432 slots vs 23,000,000 guard
GUARD_S4 = 1_470_000
GUARD_S5 = 23_000_000

def is_prime_by_trial_division(value: int) -> bool:
    """Deterministic primality check for the two ~31-bit field moduli."""
    if value < 2:
        return False
    if value % 2 == 0:
        return value == 2
    divisor = 3
    while divisor * divisor <= value:
        if value % divisor == 0:
            return False
        divisor += 2
    return True


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
            "stage": stage,
            "reason": reason,
            "processed_columns": processed,
            "total_columns": total,
            "elapsed_seconds_process_time": round(elapsed, 6),
            "peak_rss_bytes": rss,
            "peak_rss_mib": round(rss / 1024**2, 3),
            "columns_per_second": (round(processed / elapsed, 3)
                                   if elapsed > 0 else None),
            "share_of_class_not_enumerated": (round(1.0 - processed / total, 9)
                                              if total else None),
        }


# ---------------------------------------------------------------------------
# Finite geometry: the 3x3x3 box B_2 inside the extended box C_2.
# ---------------------------------------------------------------------------

BASE_SITES: tuple[tuple[int, int, int], ...] = tuple(itertools.product(range(3), repeat=3))
EXT_SITES: tuple[tuple[int, int, int], ...] = tuple(itertools.product(range(-1, 4), repeat=3))
EXT_INDEX = {site: i for i, site in enumerate(EXT_SITES)}
BASE_BITS = tuple(1 << EXT_INDEX[site] for site in BASE_SITES)
EXT_MASK = (1 << len(EXT_SITES)) - 1
NEIGHBORS: dict[int, tuple[int, ...]] = {}
for _site, _bit in zip(BASE_SITES, BASE_BITS):
    _nearby = []
    for _axis in range(3):
        for _sign in (-1, 1):
            _n = list(_site)
            _n[_axis] += _sign
            _nearby.append(1 << EXT_INDEX[tuple(_n)])
    NEIGHBORS[_bit] = tuple(_nearby)


def anchor_key(x_mask: int, z_mask: int) -> int:
    """Injective integer key of the translation orbit of the word X^x Z^z.

    The canonical representative is the unique translate whose support has
    coordinatewise minimum 0.  Every site of a class output word lies in
    C_2 = {-1..3}^3, so relative coordinates dx, dy, dz lie in 0..4 and the
    per-site digit ((dx*5+dy)*5+dz)*4 + (code-1) lies in 0..498 < 500.
    A word of support size n has n sorted digits, so its packed part is
    < 500^n <= 500^6 < 2^54, and key = (n << 54) | packed is injective.
    Ordering by key is exactly lexicographic (size, packed digits).
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
        vals.append((((sx - mnx) * 5 + (sy - mny)) * 5 + (sz - mnz)) * 4 + (code - 1))
    vals.sort()
    packed = 0
    for v in vals:
        packed = packed * 500 + v
    return (size << 54) | packed


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


def anchored_subset_count_bruteforce(smax: int) -> dict[int, int]:
    """A_k = k-subsets of the 27 sites touching all three min faces."""
    counts: dict[int, int] = {}
    for k in range(1, smax + 1):
        n = 0
        for combo in itertools.combinations(range(27), k):
            sites = [BASE_SITES[c] for c in combo]
            if (min(s[0] for s in sites) == 0
                    and min(s[1] for s in sites) == 0
                    and min(s[2] for s in sites) == 0):
                n += 1
        counts[k] = n
    return counts


def anchored_subset_count_inclusion_exclusion(k: int) -> int:
    """A_k = C(27,k) - 3 C(18,k) + 3 C(12,k) - C(8,k).

    Missing the x=0 face confines support to {1,2}x{0,1,2}^2 (18 sites);
    two specific faces to {1,2}^2 x {0,1,2} (12 sites); all three to
    {1,2}^3 (8 sites).
    """
    return (math.comb(27, k) - 3 * math.comb(18, k)
            + 3 * math.comb(12, k) - math.comb(8, k))


def rank_S_closed_form(smax: int) -> tuple[int, dict[str, list[int]]]:
    """rank S_C = 1 + sum_k A_k 3^k (distinct anchored orbits of class words).

    Inclusion-exclusion and brute enumeration must agree per size.
    """
    brute = anchored_subset_count_bruteforce(smax)
    ranks = 1
    for k in range(1, smax + 1):
        ie = anchored_subset_count_inclusion_exclusion(k)
        assert ie == brute[k], (k, ie, brute[k])
        ranks += brute[k] * 3**k
    return ranks, {
        "anchored_subsets_by_size_bruteforce": [brute[k] for k in range(1, smax + 1)],
        "anchored_subsets_by_size_inclusion_exclusion": [
            anchored_subset_count_inclusion_exclusion(k) for k in range(1, smax + 1)
        ],
    }


class LeadTable:
    """Fixed-capacity linear-probing map key -> pivot index.

    Key 0 (the empty word) is the empty sentinel; it can never be a lead
    because D maps nonempty words to nonempty words and kills the identity.
    Memory: 12 bytes per slot.  Capacity is predeclared with substantial
    headroom below the table-full failure mode.
    """

    __slots__ = ("cap", "mask", "keys", "vals")

    def __init__(self, cap_bits: int) -> None:
        self.cap = 1 << cap_bits
        self.mask = self.cap - 1
        self.keys = array("q", bytes(8 * self.cap))
        self.vals = array("i", bytes(4 * self.cap))

    def slot_of(self, key: int) -> int:
        """Slot holding key, or the empty slot where key belongs."""
        slot = (key * 0x9E3779B97F4A7C15) & 0x1FFFFFFFFFFFFF & self.mask
        keys = self.keys
        while True:
            k = keys[slot]
            if k == key or k == 0:
                return slot
            slot = (slot + 1) & self.mask


# ---------------------------------------------------------------------------
# The exact elimination (the certificate).
# ---------------------------------------------------------------------------

def eliminate(stage: str, smax: int, prime: int, budget_s: float,
              cap_bits: int, guard_pivots: int,
              log_interval: int = RSS_CHECK_COLUMN_INTERVAL) -> dict[str, object]:
    """Stream all class columns; return the exact F_p rank and measurements.

    Column echelon: the working vector w maps orbit-key -> coefficient in
    F_p.  The lead is max(w) (static (size, blob)-descending priority).
    On a fresh lead the normalized vector MINUS its lead entry is appended
    to the CSR store; on an existing lead the pivot row is subtracted from
    w in place (no copy: pivot rows never contain their own lead).  The
    number of pivots equals the F_p column rank regardless of lead order.

    Raises ResourceWall on process_time budget or RSS breach, carrying the
    exact processed/total/elapsed/RSS numbers, checked at the predeclared
    interval and at stage end.  Asserts Lemma A on every output word.
    """
    t0 = process_time()
    gc.disable()
    lead = LeadTable(cap_bits)
    offsets = array("q", [0])
    keys = array("q")
    coeffs = array("i")
    digest = hashlib.sha256()
    pcount = stored = maxlen = reductions = ops = 0
    anchored_cols = 0
    cid = -1
    try:
        for x, z in iter_class_columns(smax):
            cid += 1
            if not x and not z:
                anchored_cols += 1  # identity column: anchored, zero image
                continue
            if _is_anchored(x | z):
                anchored_cols += 1
            w: dict[int, int] = {}
            for (nx, nz), c in image_terms(x, z).items():
                assert not (nx | nz) & ~EXT_MASK, "Lemma A violated: word left C_R"
                k = anchor_key(nx, nz)
                assert (k >> 54) <= smax + 1, "Lemma A violated: support size > s+1"
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
                    digest.update(cid.to_bytes(4, "little"))
                    digest.update(lk.to_bytes(8, "little"))
                    pcount += 1
                    stored += len(w)
                    if len(w) > maxlen:
                        maxlen = len(w)
                    if pcount > guard_pivots:
                        raise ResourceWall(
                            stage,
                            f"pivot count {pcount} exceeded plan guard {guard_pivots}",
                            cid + 1, class_column_count(smax),
                            process_time() - t0, peak_rss_bytes())
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
                print(f"  [{stage}] col {cid:>10,} piv={pcount:>10,} "
                      f"nnz={stored:>12,} avg={stored / max(pcount, 1):5.2f} "
                      f"max={maxlen} red={reductions:>10,} ops={ops / 1e9:6.3f}G "
                      f"rss={rss / 2**20:7.0f}MiB t={el:7.0f}s", flush=True)
                if el > budget_s:
                    raise ResourceWall(
                        stage, f"process_time budget {budget_s}s exceeded",
                        cid + 1, class_column_count(smax), el, rss)
                if rss > RSS_WALL_BYTES:
                    raise ResourceWall(
                        stage, f"peak RSS {rss / 2**20:.1f} MiB exceeded the "
                        f"{RSS_WALL_BYTES / 10**9:.0f} GB wall",
                        cid + 1, class_column_count(smax), el, rss)
    finally:
        gc.enable()
    el = process_time() - t0
    rss = peak_rss_bytes()
    if el > budget_s:
        raise ResourceWall(stage, f"process_time budget {budget_s}s exceeded at stage end",
                           cid + 1, class_column_count(smax), el, rss)
    if rss > RSS_WALL_BYTES:
        raise ResourceWall(stage, f"peak RSS {rss / 2**20:.1f} MiB exceeded the "
                           f"{RSS_WALL_BYTES / 10**9:.0f} GB wall at stage end",
                           cid + 1, class_column_count(smax), el, rss)
    return {
        "claim_tag": "[COMPUTATION]",
        "stage": stage,
        "smax": smax,
        "prime": prime,
        "columns": cid + 1,
        "columns_closed_form": class_column_count(smax),
        "anchored_columns": anchored_cols,
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
        "rss_measurement": "process-lifetime ru_maxrss at measurement end; conservative",
    }


def _is_anchored(support: int) -> bool:
    mnx = mny = mnz = 99
    m = support
    while m:
        bit = m & -m
        m ^= bit
        sx, sy, sz = EXT_SITES[bit.bit_length() - 1]
        if sx < mnx:
            mnx = sx
        if sy < mny:
            mny = sy
        if sz < mnz:
            mnz = sz
    return mnx == 0 and mny == 0 and mnz == 0


# ---------------------------------------------------------------------------
# Greedy static row-discovery lower bound (recorded, NOT the certificate).
# ---------------------------------------------------------------------------

def greedy_lower_bound(smax: int, rule: str, budget_s: float) -> dict[str, object]:
    """Certified-but-stalling greedy rank lower bound on the natural order.

    For each streamed column, coalesce its image by translation-orbit key.
    If none of the resulting nonzero rows is an already-claimed lead, claim
    one (rule 'max' claims the largest, 'min' the smallest).  The claimed
    keys k_1..k_T and selected columns c_1..c_T form a lower-triangular
    submatrix with nonzero diagonal: by the no-hit rule c_j's image avoids
    every k_i with i < j, while c_i's image contains k_i.  Hence rank >= T.
    The count stalls far below the exact rank because most columns eventually
    hit a claimed key; it can never exceed the rank.
    """
    t0 = process_time()
    leads: set[int] = set()
    pivots = 0
    total = 0
    cid = -1
    for x, z in iter_class_columns(smax):
        cid += 1
        if not x and not z:
            continue
        total += 1
        projected: dict[int, int] = {}
        for (nx, nz), c in image_terms(x, z).items():
            k = anchor_key(nx, nz)
            projected[k] = projected.get(k, 0) + c
        ks = {k for k, c in projected.items() if c}
        if ks & leads:
            continue
        leads.add(max(ks) if rule == "max" else min(ks))
        pivots += 1
        if cid % RSS_CHECK_COLUMN_INTERVAL == 0:
            el = process_time() - t0
            rss = peak_rss_bytes()
            if el > budget_s:
                raise ResourceWall(
                    f"greedy_{rule}", f"process_time budget {budget_s}s exceeded",
                    cid + 1, class_column_count(smax), el, rss)
            if rss > RSS_WALL_BYTES:
                raise ResourceWall(
                    f"greedy_{rule}",
                    f"peak RSS {rss / 2**20:.1f} MiB exceeded the "
                    f"{RSS_WALL_BYTES / 10**9:.0f} GB wall",
                    cid + 1, class_column_count(smax), el, rss)
    el = process_time() - t0
    rss = peak_rss_bytes()
    if el > budget_s:
        raise ResourceWall(
            f"greedy_{rule}", f"process_time budget {budget_s}s exceeded at stage end",
            cid + 1, class_column_count(smax), el, rss)
    if rss > RSS_WALL_BYTES:
        raise ResourceWall(
            f"greedy_{rule}",
            f"peak RSS {rss / 2**20:.1f} MiB exceeded the "
            f"{RSS_WALL_BYTES / 10**9:.0f} GB wall at stage end",
            cid + 1, class_column_count(smax), el, rss)
    return {
        "claim_tag": "[COMPUTATION]",
        "smax": smax,
        "rule": (f"claim {'largest' if rule == 'max' else 'smallest'} image orbit key, "
                 "natural class column order"),
        "greedy_lower_bound_T": pivots,
        "columns_seen": total,
        "elapsed_seconds_process_time": round(el, 6),
        "peak_rss_bytes": rss,
        "peak_rss_mib": round(rss / 1024**2, 3),
        "certified_status": ("valid rank lower bound by lower-triangular minor; "
                             "stalls strictly below the exact rank"),
    }


# ---------------------------------------------------------------------------
# Defect traps (cheap, decisive, run before anything expensive).
# ---------------------------------------------------------------------------

def defect_traps() -> dict[str, object]:
    origin = 1 << EXT_INDEX[(0, 0, 0)]
    dz = image_terms(0, origin)
    dz_ok = dz == {(origin, origin): 1}  # D(Z_0) = +X_0 Z_0
    dx = image_terms(origin, 0)
    dx_ok = len(dx) == 6 and set(dx.values()) == {-1}  # six bond terms, all -1
    # pi D h = 0 on orbit sums, h = X_0 + sum_i Z_0 Z_{e_i} (a=b=1).
    words: dict[tuple[int, int], int] = {(origin, 0): 1}
    for axis in range(3):
        s = [0, 0, 0]
        s[axis] = 1
        nb = 1 << EXT_INDEX[tuple(s)]
        key = (0, origin | nb)
        words[key] = words.get(key, 0) + 1
    orbit_sum: dict[int, int] = {}
    for (x, z), coeff in words.items():
        for (nx, nz), c in image_terms(x, z).items():
            k = anchor_key(nx, nz)
            orbit_sum[k] = orbit_sum.get(k, 0) + coeff * c
    pdih_ok = not any(orbit_sum.values())
    return {
        "D_of_Z0_is_plus_X0Z0": dz_ok,
        "D_of_X0_has_six_bond_terms_all_minus1": dx_ok,
        "pi_D_h_zero_on_orbit_sums": pdih_ok,
        "trap_detail": {
            "D_Z0_terms": {f"{k[0]:x},{k[1]:x}": v for k, v in dz.items()},
            "D_X0_term_count": len(dx),
            "D_X0_coefficients": sorted(set(dx.values())),
        },
    }


# ---------------------------------------------------------------------------
# Orchestration.
# ---------------------------------------------------------------------------

def check(name: str, passed: bool, detail: str) -> dict[str, object]:
    return {"name": name, "passed": bool(passed), "detail": detail}


def load_prior_baseline() -> dict[str, object]:
    if not PRIOR_RESULT.exists():
        return {}
    payload = json.loads(PRIOR_RESULT.read_text(encoding="utf-8"))
    wall = payload.get("data", {}).get("observed_wall_probe", {}).get("resource_wall", {})
    return {
        "prior_front_observed_wall": {
            "processed_columns": wall.get("processed_columns"),
            "total_columns": wall.get("total_columns"),
            "elapsed_seconds": wall.get("elapsed_seconds"),
            "peak_rss_mib": wall.get("peak_rss_mib"),
            "measured_columns_per_second": wall.get("measured_columns_per_second"),
        },
        "statement": ("wave-13 row-oriented streamed-CSR baseline (stores every column "
                      "image) that this column-oriented run must beat"),
    }


def main() -> int:
    started = process_time()

    # ---- Stage 0: closed-form arithmetic, cross-checked, plus traps. ----
    t0 = process_time()
    rankS4, breakdown4 = rank_S_closed_form(4)
    rankS5, breakdown5 = rank_S_closed_form(5)
    cols4 = class_column_count(4)
    cols5 = class_column_count(5)
    traps = defect_traps()
    primes_ok = all(is_prime_by_trial_division(p) for p in PRIMES)
    traps_ok = all(v for k, v in traps.items() if isinstance(v, bool))
    arithmetic_ok = (
        cols4 == CERTIFIED_S4_COLUMNS and rankS4 == CERTIFIED_S4_RANK_S
        and cols5 == S5_COLUMNS and rankS5 == S5_RANK_S and traps_ok and primes_ok
    )
    if not arithmetic_ok or process_time() - t0 > PREFLIGHT_BUDGET_S:
        print("FAIL: preflight closed-form arithmetic or defect traps")
        return 1
    print(f"[preflight] rankS(s<=4)={rankS4:,} rankS(s<=5)={rankS5:,} "
          f"cols={cols4:,}/{cols5:,} traps ok; target rank s<=5 = {S5_TARGET_RANK:,}",
          flush=True)
    # ---- Stage 1: mandatory s<=4 regression, two primes, BEFORE the run. ----
    regression: list[dict[str, object]] = []
    reg_wall: dict[str, object] | None = None
    try:
        for i, p in enumerate(PRIMES):
            r = eliminate(f"regression_s4_p{i + 1}", 4, p,
                          REGRESSION_BUDGET_S / len(PRIMES), CAP_BITS_S4, GUARD_S4)
            assert r["columns"] == CERTIFIED_S4_COLUMNS
            assert r["anchored_columns"] == CERTIFIED_S4_RANK_S
            assert r["rank_Fp"] <= CERTIFIED_S4_RANK_M, (
                "modular rank violates the analytic {I,h} upper bound")
            regression.append(r)
            print(f"[regression s<=4 p={p}] rank_Fp={r['rank_Fp']:,} "
                  f"in {r['elapsed_seconds_process_time']}s "
                  f"rss={r['peak_rss_mib']}MiB", flush=True)
    except ResourceWall as wall_exc:
        reg_wall = wall_exc.record
    if reg_wall is not None or any(r["rank_Fp"] != CERTIFIED_S4_RANK_M for r in regression):
        write_payload(regression, [], None, None, None, None, rankS4, rankS5,
                      cols4, cols5, breakdown4, breakdown5, traps, reg_wall,
                      started, verdict="regression_failed")
        print("FAIL: s<=4 regression did not reproduce the certified values")
        return 1

    # ---- Stage 2: greedy lower bounds on the calibration class. ----
    greedy: list[dict[str, object]] = []
    greedy_wall: dict[str, object] | None = None
    try:
        for rule in ("max", "min"):
            g = greedy_lower_bound(4, rule, GREEDY_BUDGET_S / 2)
            greedy.append(g)
            print(f"[greedy s<=4 rule={rule}] T={g['greedy_lower_bound_T']:,} "
                  f"in {g['elapsed_seconds_process_time']}s", flush=True)
    except ResourceWall as wall_exc:
        greedy_wall = wall_exc.record

    # ---- Stage 3: the headline s<=5 class, two primes. ----
    headline: list[dict[str, object]] = []
    headline_wall: dict[str, object] | None = None
    try:
        for i, p in enumerate(PRIMES):
            r = eliminate(f"headline_s5_p{i + 1}", 5, p, S5_PRIME_BUDGET_S,
                          CAP_BITS_S5, GUARD_S5)
            assert r["columns"] == S5_COLUMNS
            assert r["anchored_columns"] == S5_RANK_S
            assert r["rank_Fp"] <= S5_TARGET_RANK, (
                "modular rank violates the analytic {I,h} upper bound")
            headline.append(r)
            print(f"[headline s<=5 p={p}] rank_Fp={r['rank_Fp']:,} "
                  f"in {r['elapsed_seconds_process_time']}s "
                  f"rss={r['peak_rss_mib']}MiB", flush=True)
    except ResourceWall as wall_exc:
        headline_wall = wall_exc.record

    failed = write_payload(regression, greedy, headline, greedy_wall, headline_wall,
                           load_prior_baseline(), rankS4, rankS5, cols4, cols5,
                           breakdown4, breakdown5, traps, reg_wall, started,
                           verdict=None)
    return 1 if failed else 0


def write_payload(regression, greedy, headline, greedy_wall, headline_wall,
                  baseline, rankS4, rankS5, cols4, cols5, breakdown4,
                  breakdown5, traps, reg_wall, started, verdict) -> list[str]:
    """Assemble the result JSON, the checks, and the honest verdict."""
    completed_s5 = [r for r in (headline or []) if r.get("rank_Fp") is not None]
    ranks_s5 = sorted({r["rank_Fp"] for r in completed_s5})
    quotient_certified = (len(completed_s5) == len(PRIMES)
                          and ranks_s5 == [S5_TARGET_RANK])
    regression_ok = (reg_wall is None and len(regression) == len(PRIMES)
                     and all(r["rank_Fp"] == CERTIFIED_S4_RANK_M for r in regression))

    if headline_wall is not None or len(completed_s5) < len(PRIMES):
        outcome = "observed_wall" if headline_wall is not None else "incomplete"
        quotient = None
    elif ranks_s5 == [S5_TARGET_RANK]:
        outcome = "certified_quotient_1"
        quotient = 1
    elif len(ranks_s5) == 1:
        outcome = "modular_rank_below_analytic_bound"
        quotient = None
    else:
        outcome = "prime_disagreement"
        quotient = None
    if verdict == "regression_failed":
        outcome = "regression_failed"

    budgets_respected = all(
        r["elapsed_seconds_process_time"] <= r["budget_seconds_process_time"]
        and r["peak_rss_bytes"] <= RSS_WALL_BYTES
        for r in list(regression) + list(completed_s5)
    )

    greedy_records = list(greedy or [])
    greedy_best = max((g["greedy_lower_bound_T"] for g in greedy_records), default=None)
    greedy_ratio = (round(greedy_best / CERTIFIED_S4_RANK_M, 6)
                    if greedy_best is not None else None)

    checks = [
        check("closed_form_arithmetic_and_traps",
              cols4 == CERTIFIED_S4_COLUMNS and rankS4 == CERTIFIED_S4_RANK_S
              and cols5 == S5_COLUMNS and rankS5 == S5_RANK_S
              and all(v for k, v in traps.items() if isinstance(v, bool)),
              "column counts and rank S closed forms reproduce the frozen values; "
              "inclusion-exclusion equals brute anchored enumeration; D/Z0, D/X0 and "
              "pi D h traps all pass"),
        check("regression_s4_both_primes_reproduced", regression_ok,
              "the s<=4 class reproduces columns 1,503,766, anchored census 822,334 "
              "and rank 822,332 at both primes before the new class ran"),
        check("budgets_and_rss_respected_on_completed_stages", budgets_respected,
              "every completed stage stayed inside its predeclared process_time "
              "budget and the 6 GB RSS wall"),
        check("greedy_gap_recorded",
              greedy_best is not None and greedy_best < CERTIFIED_S4_RANK_M,
              f"greedy static lower bound stalls at T={greedy_best} "
              f"(ratio {greedy_ratio}) of the certified exact rank on the s<=4 "
              "calibration class; the exact elimination is the certificate"),
        check("headline_outcome_recorded_honestly",
              outcome in {"certified_quotient_1", "observed_wall",
                          "modular_rank_below_analytic_bound",
                          "regression_failed", "incomplete", "prime_disagreement"}
              and (quotient == 1) == quotient_certified,
              f"s<=5 outcome is {outcome}; a quotient is claimed only when both "
              "primes completed at the analytic bound"),
        check("full_box_z3_r2_remains_undecided_statement", True,
              "the unresolved section restates that the FULL Z^3 radius-2 box is "
              "NOT decided; only the support-size class s<=5 is at stake here"),
    ]

    headline_block: dict[str, object] = {
        "class": "C(3,2,5): rational span of ordered-Pauli words supported in "
                 "B_2={0,1,2}^3 with word support size <= 5",
        "columns_closed_form": S5_COLUMNS,
        "rank_S_closed_form": S5_RANK_S,
        "analytic_upper_bound_rank_M": S5_TARGET_RANK,
        "target_quotient_dimension": 1,
        "outcome": outcome,
        "cases": completed_s5,
    }
    if headline_wall is not None:
        headline_block["observed_wall"] = headline_wall
        headline_block["statement"] = (
            "The s<=5 elimination was stopped by an OBSERVED budget breach; no rank "
            "or quotient is claimed. The wall record carries the exact share not "
            "enumerated.")
    elif outcome == "certified_quotient_1":
        headline_block["certified_quotient_dimension"] = 1
        headline_block["certification"] = (
            "rank_Fp = rank S - 2 at both primes with rank_Fp <= rank_Q <= rank S - 2 "
            "forces rank_Q = 14,757,410; quotient = rank S - rank M - 1 = 1: modulo "
            "finite translation shifts and the identity, the only class density with "
            "vanishing projected commutator is the Hamiltonian density h")
        headline_block["kernel_nullity"] = S5_COLUMNS - S5_TARGET_RANK
        headline_block["trivial_divergence_dimension"] = S5_COLUMNS - S5_RANK_S
    elif outcome == "modular_rank_below_analytic_bound":
        headline_block["statement"] = (
            f"both primes completed and agree at rank {ranks_s5[0]:,} < "
            f"{S5_TARGET_RANK:,}; this is only a lower bound on rank_Q and no "
            "quotient is certified")
    elif outcome == "prime_disagreement":
        headline_block["statement"] = "the two primes disagree; no claim is made"

    greedy_block: dict[str, object] = {
        "claim_tag": "[COMPUTATION]",
        "calibration_class": "C(3,2,4)",
        "exact_rank": CERTIFIED_S4_RANK_M,
        "records": greedy_records,
        "best_greedy_lower_bound": greedy_best,
        "best_greedy_ratio_of_exact": greedy_ratio,
        "statement": (
            f"the certified greedy lower bound stalls at {greedy_ratio} of the exact "
            "rank on the calibration class, so the greedy certificate cannot reach "
            "rank S - 2 and the exact column echelon is the certificate"),
    }
    if greedy_wall is not None:
        greedy_block["observed_wall"] = greedy_wall

    unresolved: list[dict[str, object]] = [
        {
            "claim_tag": "[UNRESOLVED]",
            "scope": "full Z^3 radius-2 box",
            "statement": (
                "The FULL Z^3 radius-2 (3x3x3) ordered-Pauli box is NOT decided by "
                "this artifact and was not decided before it: only the support-size "
                "class s<=5 is at stake. Densities using any word with six or more "
                "nonidentity sites in that box remain unchecked."),
        },
        {
            "claim_tag": "[UNRESOLVED]",
            "scope": "classes s>=6 and all larger boxes",
            "statement": ("support-size classes s>=6 in Z^3 R=2, the full boxes "
                          "Z^3 R>=2 and Z^2 R>=3, quasilocal and non-translation-"
                          "covariant charges remain undecided"),
        },
        {
            "claim_tag": "[UNRESOLVED]",
            "scope": "coupling ratio",
            "statement": ("all certified values are a=b=1 statements; exceptional "
                          "nonzero coupling ratios are not excluded"),
        },
        {
            "claim_tag": "[UNRESOLVED]",
            "scope": "integrability",
            "statement": ("nothing here proves non-integrability, a solution, or an "
                          "all-size theorem for the three-dimensional Ising model"),
        },
    ]
    if outcome in {"observed_wall", "incomplete"}:
        unresolved.insert(0, {
            "claim_tag": "[UNRESOLVED]",
            "scope": "Z^3 R=2 class s<=5",
            "statement": (
                "the s<=5 elimination hit an OBSERVED resource wall; the wall record "
                "states the exact processed/total share, process_time and peak RSS; "
                "no rank or quotient is claimed for this class"),
        })
    elif outcome == "modular_rank_below_analytic_bound":
        unresolved.insert(0, {
            "claim_tag": "[UNRESOLVED]",
            "scope": "Z^3 R=2 class s<=5",
            "statement": (f"modular rank {ranks_s5[0]:,} < rank S - 2; only the "
                          "two-sided bound below is certified"),
        })

    payload = {
        "provenance": {
            "script": SCRIPT,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": ".venv/bin/python",
            "python": platform.python_version(),
            "platform": platform.platform(),
            "method": (
                "column-oriented exact elimination over F_p with streamed class "
                "columns, compact CSR pivot store excluding lead entries, "
                "open-addressing flat-array lead map, injective size<<54|base-500 "
                "orbit keys, static descending (size, blob) lead priority; "
                "certificate = wave-13 sandwich against the analytic {I,h} kernel"),
            "certifying_arithmetic": ("Python integers; finite-field elimination "
                                      "modulo the listed primes; no floats"),
            "clock": "time.process_time() only; never wall clock or signal alarms",
            "budgets_seconds_process_time": {
                "preflight": PREFLIGHT_BUDGET_S,
                "regression_total": REGRESSION_BUDGET_S,
                "greedy_total": GREEDY_BUDGET_S,
                "s5_per_prime": S5_PRIME_BUDGET_S,
            },
            "rss_wall_bytes": RSS_WALL_BYTES,
            "rss_measurement": "process-lifetime ru_maxrss at measurement end; conservative",
            "baseline_beaten": baseline,
            "total_elapsed_seconds_process_time": round(process_time() - started, 6),
            "peak_rss_bytes": peak_rss_bytes(),
            "peak_rss_mib": round(peak_rss_bytes() / 1024**2, 3),
        },
        "data": {
            "claim_tags": ["[THEOREM]", "[LEMMA]", "[COMPUTATION]", "[EXTERNAL]", "[UNRESOLVED]"],
            "construction": {
                "claim_tag": "[LEMMA]",
                "ordered_pauli_convention": "X^a Z^b",
                "generator_orientation": "D = (1/2)[H, .] with H = sum X + sum Z Z over nearest neighbours, a=b=1",
                "ad_formulas": (
                    "(1/2)[X_s,X^aZ^b]=[b_s=1]X^(a+e_s)Z^b; "
                    "(1/2)[Z_uZ_v,X^aZ^b]=-[a_u xor a_v=1]X^aZ^(b+e_u+e_v)"),
                "class": "C(3,2,s) = span of words supported in B_2 with word support size <= s",
                "class_quotient": (
                    "dim ker(pi o D|_C)/(C cap ker pi + Q I) = rank S_C - rank M_C - 1 "
                    "(Lemmas A/B of proofs/extensive_3d_r2.md)"),
                "sandwich": "rank_Fp(M_C) <= rank_Q(M_C) <= rank S_C - 2 (kernel I and h)",
                "orbit_key_packing": (
                    "key = size<<54 | packed, digits ((dx*5+dy)*5+dz)*4+code-1 < 500 "
                    "sorted, at most 6 digits; injective since packed < 500^6 < 2^54"),
            },
            "closed_form_arithmetic": {
                "s4": {"columns": cols4, "rank_S": rankS4, **breakdown4},
                "s5": {"columns": cols5, "rank_S": rankS5, **breakdown5},
                "defect_traps": traps,
            },
            "regression_s4": {
                "certified_values": {
                    "columns": CERTIFIED_S4_COLUMNS,
                    "rank_S": CERTIFIED_S4_RANK_S,
                    "rank_M_both_primes": CERTIFIED_S4_RANK_M,
                },
                "cases": regression,
                "reproduced": regression_ok,
                "observed_wall": reg_wall,
                "statement": ("mandatory preflight: the certified s<=4 values were "
                              "reproduced BEFORE the s<=5 class ran"),
            },
            "greedy_gap": greedy_block,
            "headline_s5": headline_block,
            "unresolved": unresolved,
        },
        "checks": checks,
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    failed = [c["name"] for c in checks if not c["passed"]]
    if failed:
        print("FAIL: " + ", ".join(failed))
    else:
        print("checks: all passed; s<=5 outcome = " + outcome)
    return failed


if __name__ == "__main__":
    raise SystemExit(main())
