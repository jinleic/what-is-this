#!/usr/bin/env python3
"""Exact merged height-schedule sieve union at n=36 (wave 14).

[THEOREM] This producer recomputes the complete finite union of the 41 wave-13
height-schedule families together with all 363 admissible two-block schedules with
heights prescribed at k>=2 cut positions, unions them with the certified
L, C, F2, Q2, F3, Q3 defect families, and evaluates the union by an exact
Mobius sieve over the closure of admissible cut sequences.  Every finite count is a
Python integer.  The only rounding is the directed 40-decimal-place floor of the
120-dps interval enclosure of the resulting Ising lower endpoint.

The producer replays all decisive computation with stakes:
  * the raw positive-height profile is re-derived by a packed visited-set DFS;
  * the admissible rise table for block lengths 5..13 is enumerated and every one
    of the 392 admissible two-block schedules is recorded (not just winners);
  * the height-parity rule q == l (mod 2) is verified for every entry;
  * all 820 pairwise old-grid separators are re-derived;
  * all 81406 member pairs are checked for consistent (intersecting) unions;
  * the sieve coefficient of every active term is computed once by the hitting-mask
    rule and independently by the explicit subset-union Mobius sum;
  * every active term's N, C, L, C\\cap L, F2, Q2, F3, Q3 intersections
    are computed by exact subset-NFA or product DPs and stored row by row;
  * the old-41 projection reproduces the certified wave-13 union exactly.

Run:
    PYTHONPATH=src .venv/bin/python experiments/e141_saw_union4.py
    PYTHONPATH=src .venv/bin/python tests/test_saw_union4.py
"""
from __future__ import annotations

from bisect import bisect_right
from collections import Counter, defaultdict
from datetime import datetime, timezone
from decimal import Decimal, ROUND_FLOOR, localcontext
from itertools import product
import hashlib
import json
from pathlib import Path
import resource
import sys
import time

import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = "experiments/e141_saw_union4.py"
RESULT = ROOT / "results" / "bounds" / "saw_union4.json"
SAW_SOURCE = ROOT / "sources" / "fulltext" / "schram_barkema_bisseling2011.pdf"

C36 = 2941370856334701726560670
C36_SHA256 = "898d8333e8d70acd279e29f226ff108550b132463aa662f86df590b0ad72cb12"
WAVE13_UNION = 507841195880595165555
DPS = 120
REPORT_PLACES = 40
R = 35
M = (1, 2, 4)  # {-e1, +e2, +e3}; all raise h.
STEPS = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))
DH = (-1, 1, 1, -1, 1, -1)
OFF = 40

CANONICAL_TAIL = [1, 3, 9, 45, 171, 837, 3411, 16425, 69525, 331947,
                  1436643, 6827355, 29971017, 142016661]
CANONICAL_BLOCKS = {
    1: {1: 3}, 2: {2: 9}, 3: {3: 27}, 4: {4: 81},
    5: {3: 108, 5: 243}, 6: {4: 648, 6: 729},
    7: {3: 432, 5: 2916, 7: 2187}, 8: {4: 5778, 6: 11664, 8: 6561},
    9: {3: 1566, 5: 34668, 7: 43740, 9: 19683},
    10: {4: 48168, 6: 167670, 8: 157464, 10: 59049},
    11: {3: 5832, 5: 411156, 7: 729000, 9: 551124, 11: 177147},
    12: {4: 394524, 6: 2356128, 8: 2969946, 10: 1889568, 12: 531441},
    13: {3: 20952, 5: 4815450, 7: 11544444, 9: 11573604, 11: 6377292,
         13: 1594323},
}

checks: list[dict[str, object]] = []


def check(name: str, passed: bool, detail: str) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})
    if not passed:
        print(f"FAIL {name}: {detail}")
        raise AssertionError(name)
    print(f"PASS {name}: {detail}")


def peak_rss_bytes() -> int:
    raw = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(raw if sys.platform == "darwin" else raw * 1024)


# ---------------------------------------------------------------------------
# Raw positive-height profile by packed-coordinate visited-set DFS.
# ---------------------------------------------------------------------------


def pack(x: int, y: int, z: int) -> int:
    return ((x + OFF) << 14) | ((y + OFF) << 7) | (z + OFF)


def pack_delta(dx: int, dy: int, dz: int) -> int:
    return (dx << 14) + (dy << 7) + dz


PACKED_STEPS = tuple(pack_delta(*step) for step in STEPS)


def positive_height_profile(max_depth: int) -> tuple[list[int], list[dict[int, int]]]:
    """Tail T(n) and strict-record block counts B(n,q) through ``max_depth``."""
    tail = [0] * (max_depth + 1)
    blocks = [defaultdict(int) for _ in range(max_depth + 1)]
    origin = pack(0, 0, 0)
    visited = {origin}

    def visit(point: int, height: int, previous_max: int, depth: int) -> None:
        tail[depth] += 1
        if depth and height > previous_max:
            blocks[depth][height] += 1
        if depth == max_depth:
            return
        next_previous_max = max(previous_max, height)
        for letter, delta in enumerate(PACKED_STEPS):
            nxt = point + delta
            if nxt in visited:
                continue
            next_height = height + DH[letter]
            if next_height <= 0:
                continue
            visited.add(nxt)
            visit(nxt, next_height, next_previous_max, depth + 1)
            visited.remove(nxt)

    visit(origin, 0, 0, 0)
    return tail, [dict(row) for row in blocks]


# ---------------------------------------------------------------------------
# Height schedules, cut sequences, and the member families.
# ---------------------------------------------------------------------------


def sched2cuts(s: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    """Absolute cut positions and heights of a block-length schedule."""
    if len(s) == 5:
        l1, q1, l2, q2, _t = s
        return ((l1, q1), (l1 + l2, q1 + q2))
    l1, q1, l2, q2, l3, q3, _t = s
    return ((l1, q1), (l1 + l2, q1 + q2), (l1 + l2 + l3, q1 + q2 + q3))


def member_count(cuts: tuple[tuple[int, int], ...]) -> int:
    """|F_cuts|: concatenated strict-record blocks and a positive tail."""
    prevp, prevh, total = 0, 0, 1
    for p, h in cuts:
        block = CANONICAL_BLOCKS.get(p - prevp, {}).get(h - prevh, 0)
        total *= block
        if total == 0:
            return 0
        prevp, prevh = p, h
    return total * CANONICAL_TAIL[R - prevp]


def allowed_heights_sched(schedule: tuple[int, ...], position: int) -> set[int]:
    """Exact finite allowed-height set at a step, bounded by R."""
    if len(schedule) == 5:
        l1, q1, l2, q2, _t = schedule
        if position < l1:
            return set(range(1, q1))
        if position == l1:
            return {q1}
        if position < l1 + l2:
            return set(range(q1 + 1, q1 + q2))
        if position == l1 + l2:
            return {q1 + q2}
        return set(range(q1 + q2 + 1, R + 1))
    l1, q1, l2, q2, l3, q3, _t = schedule
    cut2, cut3 = l1 + l2, l1 + l2 + l3
    if position < l1:
        return set(range(1, q1))
    if position == l1:
        return {q1}
    if position < cut2:
        return set(range(q1 + 1, q1 + q2))
    if position == cut2:
        return {q1 + q2}
    if position < cut3:
        return set(range(q1 + q2 + 1, q1 + q2 + q3))
    if position == cut3:
        return {q1 + q2 + q3}
    return set(range(q1 + q2 + q3 + 1, R + 1))


def first_separator(left: tuple[int, ...], right: tuple[int, ...]) -> dict[str, object] | None:
    """[LEMMA] an empty allowed-height slice proves H_left and H_right disjoint."""
    for position in range(1, R + 1):
        lh = allowed_heights_sched(left, position)
        rh = allowed_heights_sched(right, position)
        if not lh.intersection(rh):
            return {"step": position, "left_allowed": sorted(lh), "right_allowed": sorted(rh)}
    return None


def build_members(blocks: dict[int, dict[int, int]]) -> tuple[dict[tuple[tuple[int, int], ...], str],
                                                              list[tuple[str, tuple[int, ...]]]]:
    """All 404 members: the wave-13 41 plus the admissible two-block sweep.

    The two-block sweep enumerates the full table l_1, l_2, t in [5..13]
    with t = 35 - l_1 - l_2 and all admissible rises; 392 schedules are
    admissible of which 29 already occur in the wave-13 grid, giving
    41 + 392 - 29 = 404 members.
    """
    members: dict[tuple[tuple[int, int], ...], str] = {}
    for q1 in (3, 5, 7, 9, 11):
        for q2 in (3, 5, 7, 9, 11):
            members.setdefault(sched2cuts((11, q1, 11, q2, 13)), f"H{q1}_{q2}")
    for q in (5, 7, 9, 11):
        members.setdefault(sched2cuts((13, q, 11, 3, 11)), f"V{q}")
    for q1 in (3, 5, 7):
        for q3 in (3, 5, 7, 9):
            members.setdefault(sched2cuts((7, q1, 7, 3, 9, q3, 12)), f"J{q1}_{q3}")
    twoblock: list[tuple[int, ...]] = []
    for l1 in range(5, 14):
        for l2 in range(5, 14):
            t = R - l1 - l2
            if not 5 <= t <= 13:
                continue
            for q1 in blocks[l1]:
                for q2 in blocks[l2]:
                    twoblock.append((l1, q1, l2, q2, t))
                    members.setdefault(sched2cuts((l1, q1, l2, q2, t)),
                                   f"N{l1}_{q1}_{l2}_{q2}_{t}")
    return members, [(f"N{s[0]}_{s[1]}_{s[2]}_{s[3]}_{s[4]}", s) for s in twoblock]


def schedule_specs_of(members: dict[tuple[tuple[int, int], ...], str]) -> list[tuple[str, tuple[int, ...]]]:
    specs: list[tuple[str, tuple[int, ...]]] = []
    for q1 in (3, 5, 7, 9, 11):
        for q2 in (3, 5, 7, 9, 11):
            specs.append((f"H{q1}_{q2}", (11, q1, 11, q2, 13)))
    for q in (5, 7, 9, 11):
        specs.append((f"V{q}", (13, q, 11, 3, 11)))
    for q1 in (3, 5, 7):
        for q3 in (3, 5, 7, 9):
            specs.append((f"J{q1}_{q3}", (7, q1, 7, 3, 9, q3, 12)))
    return specs


# ---------------------------------------------------------------------------
# Closure of admissible cut sequences and the exact sieve.
# ---------------------------------------------------------------------------


def closure_and_sieve(members: dict[tuple[tuple[int, int], ...], str],
                     blocks: dict[int, dict[int, int]]
                     ) -> tuple[list[tuple[tuple[int, int], ...]],
                              dict[tuple[tuple[int, int], ...], tuple[int, list[str]]],
                              dict[int, int], int]:
    POS = sorted({p for cuts in members for p, _h in cuts})
    member_sets = [frozenset(m) for m in members]
    member_names = list(members.values())
    cands: list[tuple[tuple[int, int], ...]] = []

    def transitions(p: int, h: int):
        for pp in POS:
            if pp <= p:
                continue
            dp = pp - p
            if dp > 13:
                break
            for dh in blocks.get(dp, {}):
                yield pp, h + dh

    def extend(cuts: tuple[tuple[int, int], ...], p: int, h: int) -> None:
        for pp, hp in transitions(p, h):
            new = cuts + ((pp, hp),)
            cands.append(new)
            extend(new, pp, hp)

    for p, q in [(p, q) for p in POS if p <= 13 for q in blocks[p]]:
        cands.append(((p, q),))
        extend(((p, q),), p, q)

    gsize = Counter()
    active: dict[tuple[tuple[int, int], ...], tuple[int, list[str]]] = {}
    for cuts in cands:
        U = frozenset(cuts)
        G = [ms for ms in member_sets if ms <= U]
        gsize[len(G)] += 1
        if not G:
            continue
        covered = frozenset()
        for ms in G:
            covered |= ms
        if covered != U:
            continue
        # hitting-mask sieve coefficient: sum over subsets T of cuts hitting
        # every member inside U of (-1)^{|T|+1}.
        ulist = sorted(U)
        mbits = [sum(1 << ulist.index(e) for e in ms) for ms in G]
        c = 0
        for mask in range(1, 1 << len(ulist)):
            if all(mb & mask for mb in mbits):
                c += 1 if (bin(mask).count("1") & 1) else -1
        if c:
            active[cuts] = (c, [member_names[i] for i, ms in enumerate(member_sets) if ms <= U])
    covered_chains = 0
    for cuts in cands:
        U = frozenset(cuts)
        G = [ms for ms in member_sets if ms <= U]
        if not G:
            continue
        if frozenset().union(*G) == U:
            covered_chains += 1
    return cands, active, dict(sorted(gsize.items())), covered_chains


def explicit_subset_coefficient(U: frozenset[tuple[int, int]],
                             member_sets: list[frozenset[tuple[int, int]]]) -> int:
    """Independent coefficient: sum over S subset G(U) with union S == U of
    (-1)^{|S|+1}.  (Different shape from the hitting-mask rule.)"""
    G = [ms for ms in member_sets if ms <= U]
    lidx = {e: i for i, e in enumerate(sorted(U))}
    ubit = (1 << len(lidx)) - 1
    bits = [sum(1 << lidx[e] for e in m) for m in G]
    suffix = [0] * (len(bits) + 1)
    for i in range(len(bits) - 1, -1, -1):
        suffix[i] = suffix[i + 1] | bits[i]
    total = 0

    def rec(i: int, mask: int, parity: int) -> None:
        nonlocal total
        if mask == ubit:
            total += 1 if parity == 1 else -1  # sign (-1)^{|S|+1}
        if i == len(bits):
            return
        if (mask | suffix[i]) != ubit:
            return
        rec(i + 1, mask, parity)
        rec(i + 1, mask | bits[i], parity ^ 1)

    rec(0, 0, 0)
    return total


# ---------------------------------------------------------------------------
# Subset-NFA machinery for the L, C, F2, Q2, F3, Q3 intersections.
# ---------------------------------------------------------------------------


class Sched:
    __slots__ = ("cuts", "cpos", "chgt")

    def __init__(self, cuts: tuple[tuple[int, int], ...]):
        self.cuts = tuple(cuts)
        self.cpos = tuple(p for p, _h in cuts)
        self.chgt = tuple(h for _p, h in cuts)

    def allows(self, position: int, height: int) -> bool:
        if height <= 0 or height > R:
            return False
        i = bisect_right(self.cpos, position)
        if i and self.cpos[i - 1] == position:
            return height == self.chgt[i - 1]
        lo = self.chgt[i - 1] if i else 0
        hi = self.chgt[i] if i < len(self.chgt) else R + 1
        return lo < height < hi


class Nfa:
    """Determinized subset NFA; transition table built once per process."""

    def __init__(self, letters: tuple[int, ...], moves: dict[tuple[object, int], tuple[object, ...]],
                 start, accept):
        self.letters = letters
        self.start = frozenset(start)
        self.accept = frozenset(accept)
        self.trans: dict[tuple[frozenset[object], int], frozenset[object]] = {}
        seen = {self.start}
        stack = [self.start]
        while stack:
            cur = stack.pop()
            for letter in letters:
                out = set()
                for st in cur:
                    out.update(moves.get((st, letter), ()))
                nxt = frozenset(out)
                self.trans[(cur, letter)] = nxt
                if nxt and nxt not in seen:
                    seen.add(nxt)
                    stack.append(nxt)


def nfa_schedule_dp(sch: Sched, nfa: Nfa) -> int:
    dp: dict[tuple[int, frozenset[object]], int] = {(0, nfa.start): 1}
    for position in range(1, R + 1):
        nd: defaultdict[tuple[int, frozenset[object]], int] = defaultdict(int)
        for (h, sub), cnt in dp.items():
            for l in nfa.letters:
                nsub = nfa.trans[(sub, l)]
                if not nsub:
                    continue
                nh = h + DH[l]
                if not sch.allows(position, nh):
                    continue
                nd[(nh, nsub)] += cnt
        dp = nd
    return sum(cnt for (h, sub), cnt in dp.items() if sub & nfa.accept)


def corner_moves(pred: int, defect: int, forb: int | None) -> dict[tuple[int, int], tuple[int, ...]]:
    mv: dict[tuple[int, int], tuple[int, ...]] = {}
    for l in M:
        if l == pred:
            mv[(0, l)] = (1,)
            mv[(1, l)] = (1,)
        else:
            mv[(0, l)] = (0,)
            mv[(1, l)] = (0,)
        if l != forb:
            mv[(2, l)] = (3,)
        mv[(3, l)] = (3,)
    mv[(1, defect)] = (2,)
    return mv


def fold_moves(mid: int) -> dict[tuple[int, int], tuple[int, ...]]:
    mv: dict[tuple[int, int], tuple[int, ...]] = {}
    for l in M:
        if l == 1:
            mv[(0, l)] = (1,)
            mv[(1, l)] = (1,)
        elif l == mid:
            mv[(0, l)] = (0,)
            mv[(1, l)] = (0, 2)
        else:
            mv[(0, l)] = (0,)
            mv[(1, l)] = (0,)
        if l in (2, 4):
            mv[(3, l)] = (4,)
        mv[(4, l)] = (4,)
    mv[(2, 0)] = (3,)
    return mv


def cm_corner_nfa(pred: int, defect: int, cm_letters: tuple[int, ...]):
    mv = corner_moves(pred, defect, None)
    letters = tuple(sorted(set(cm_letters) | {defect}))
    trimmed = {(st, l): mv[(st, l)] for l in letters for st in range(4) if (st, l) in mv}
    return letters, trimmed, (0,), (2, 3)


def gadget_moves(kind: str, k: int):
    if kind == "fold":
        gadget, safe = (1, 2, 0), (2, 4)
    else:
        gadget, safe = (2, 4, 3), (1, 4)
    mv: dict[tuple[object, int], tuple[object, ...]] = {}
    for g in range(k + 1):
        for gp in (0, 1, 2):
            for ns in (False, True):
                st = (g, gp, ns)
                for l in range(6):
                    outs = []
                    if gp > 0:
                        if l == gadget[gp]:
                            if gp == 2:
                                if g + 1 <= k:
                                    outs.append((g + 1, 0, True))
                            else:
                                outs.append((g, gp + 1, False))
                    elif not ns or l in safe:
                        if l in M:
                            outs.append((g, 0, False))
                        if g < k and l == gadget[0]:
                            outs.append((g, 1, False))
                    if outs:
                        mv[(st, l)] = tuple(outs)
    letters = tuple(sorted(set(M) | set(gadget)))
    return letters, mv, ((0, 0, False),), ((k, 0, False), (k, 0, True))


L_NFAS: dict[str, Nfa] = {
    "mono": Nfa(M, {(0, l): (0,) for l in M}, (0,), (0,)),
    "G1": Nfa((1, 2, 3, 4), corner_moves(4, 3, 2), (0,), (2, 3)),
    "G2": Nfa((1, 2, 4, 5), corner_moves(2, 5, 4), (0,), (2, 3)),
    "G3": Nfa((1, 2, 3, 4), corner_moves(1, 3, 2), (0,), (2, 3)),
    "G4": Nfa((1, 2, 4, 5), corner_moves(1, 5, 4), (0,), (2, 3)),
    "D1": Nfa((0, 1, 2, 4), fold_moves(2), (0,), (3, 4)),
    "D2": Nfa((0, 1, 2, 4), fold_moves(4), (0,), (3, 4)),
}
CL_NFAS: dict[str, Nfa] = {
    name: Nfa(*cm_corner_nfa(pred, defect, cm)) for name, pred, defect, cm in (
        ("G1", 4, 3, (1, 4)), ("G2", 2, 5, (1, 2)),
        ("G3", 1, 3, (1, 4)), ("G4", 1, 5, (1, 2)))
}
G_NFAS: dict[tuple[str, int], Nfa] = {
    (kind, k): Nfa(*gadget_moves(kind, k)) for kind in ("fold", "corner") for k in (2, 3)
}

_block_cache: dict[tuple[int, int, tuple[int, ...]], int] = {}
_tail_cache: dict[tuple[int, tuple[int, ...]], int] = {}


def block_word_dp(length: int, rise: int, letters: tuple[int, ...]) -> int:
    key = (length, rise, letters)
    got = _block_cache.get(key)
    if got is None:
        dp = {0: 1}
        for position in range(1, length + 1):
            nd: defaultdict[int, int] = defaultdict(int)
            for h, cnt in dp.items():
                for l in letters:
                    nh = h + DH[l]
                    if position == length:
                        if nh != rise:
                            continue
                    elif not 0 < nh < rise:
                        continue
                    nd[nh] += cnt
            dp = nd
        got = dp.get(rise, 0)
        _block_cache[key] = got
    return got


def tail_word_dp(length: int, letters: tuple[int, ...]) -> int:
    key = (length, letters)
    got = _tail_cache.get(key)
    if got is None:
        dp = {0: 1}
        for _position in range(length):
            nd: defaultdict[int, int] = defaultdict(int)
            for h, cnt in dp.items():
                for l in letters:
                    nh = h + DH[l]
                    if nh > 0:
                        nd[nh] += cnt
            dp = nd
        got = sum(dp.values())
        _tail_cache[key] = got
    return got


def cm_overlap(cuts: tuple[tuple[int, int], ...]) -> int:
    """C cap F_cuts by signed-orthant Mobius product DPs."""
    comps: list[tuple[int, int]] = []
    pp, ph = 0, 0
    for p, h in cuts:
        comps.append((p - pp, h - ph))
        pp, ph = p, h
    t = R - pp
    total = 0
    for code in product((None, 0, 1), repeat=3):
        letters = tuple(axis * 2 + (0 if sign == 0 else 1)
                      for axis, sign in enumerate(code) if sign is not None)
        if not letters:
            continue
        e = sum(sign is None for sign in code)
        value = tail_word_dp(t, letters)
        for l, q in comps:
            value *= block_word_dp(l, q, letters)
            if value == 0:
                break
        total += (-1) ** e * value
    return total


def d_row(cuts: tuple[tuple[int, int], ...], sch: Sched) -> dict[str, int]:
    lrows = {name: nfa_schedule_dp(sch, nfa) for name, nfa in L_NFAS.items()}
    clv = lrows["mono"] + sum(nfa_schedule_dp(sch, nfa) for nfa in CL_NFAS.values())
    return {
        "H": member_count(cuts),
        "C": cm_overlap(cuts),
        "L": sum(lrows.values()),
        "CL": clv,
        "F2": nfa_schedule_dp(sch, G_NFAS[("fold", 2)]),
        "Q2": nfa_schedule_dp(sch, G_NFAS[("corner", 2)]),
        "F3": nfa_schedule_dp(sch, G_NFAS[("fold", 3)]),
        "Q3": nfa_schedule_dp(sch, G_NFAS[("corner", 3)]),
    }


def row_net(row: dict[str, int]) -> int:
    return row["H"] - row["C"] - row["L"] + row["CL"] - row["F2"] - row["Q2"] - row["F3"] - row["Q3"]


# ---------------------------------------------------------------------------
# Inclusion--exclusion aggregates and endpoint.
# ---------------------------------------------------------------------------


def atanh_root_interval(magnitude: int) -> tuple[str, str]:
    previous = mp.iv.dps
    try:
        mp.iv.dps = DPS
        value = mp.iv.mpf(magnitude)
        v = mp.iv.exp(mp.iv.log(value) * (mp.iv.mpf(-1) / 36))
        text = str(mp.iv.log((1 + v) / (1 - v)) / 2).strip()
        lower, upper = text[1:-1].split(",", maxsplit=1)
        return lower.strip(), upper.strip()
    finally:
        mp.iv.dps = previous


def floor40(value: str) -> str:
    with localcontext() as context:
        context.prec = 220
        return format(Decimal(value).quantize(Decimal(1).scaleb(-REPORT_PLACES), rounding=ROUND_FLOOR), "f")


def main() -> None:
    started = time.process_time()
    stages: dict[str, float] = {}

    source_hash = hashlib.sha256(SAW_SOURCE.read_bytes()).hexdigest()
    check("external_c36_source_hash", source_hash == C36_SHA256,
          "[EXTERNAL] cached Schram--Barkema--Bisseling Table-I PDF hash")

    t0 = time.process_time()
    tail, blocks = positive_height_profile(13)
    stages["positive_height_profile"] = time.process_time() - t0
    check("positive_height_profile",
          tail == CANONICAL_TAIL and all(blocks[n] == CANONICAL_BLOCKS[n] for n in range(1, 14)),
          "[COMPUTATION] packed raw DFS reproduces every stored T(n), B(n,q)")
    parity_ok = all(q % 2 == (n % 2) for n in range(1, 14) for q in blocks[n])
    check("height_parity_rule", parity_ok,
          "[LEMMA] admissible rise q at block length l has q == l (mod 2) for every table entry")

    blocks_dict = {n: blocks[n] for n in range(1, 14)}
    t0 = time.process_time()
    members, twoblock_specs = build_members(blocks_dict)
    old_specs = schedule_specs_of(members)
    check("member_table", len(members) == 404 and len(old_specs) == 41 and len(twoblock_specs) == 392,
          "[COMPUTATION] 41 wave-13 members and the complete two-block admissible sweep to 404")
    separator_rows: list[dict[str, object]] = []
    for left_index, (left_id, left_schedule) in enumerate(old_specs):
        for right_id, right_schedule in old_specs[left_index + 1:]:
            separator = first_separator(left_schedule, right_schedule)
            if separator is None:
                raise AssertionError(f"unseparated old schedules {left_id},{right_id}")
            separator_rows.append({"left": left_id, "right": right_id, **separator})
    check("old41_pair_separators", len(separator_rows) == 820,
          "[LEMMA] every one of the 820 wave-13 schedule pairs has an exact disjoint height slice")
    stages["members_and_separators"] = time.process_time() - t0

    t0 = time.process_time()
    cands, active, gsize_census, covered_chains = closure_and_sieve(members, blocks_dict)
    coeff_dist = dict(sorted(Counter(c for c, _ in active.values()).items()))
    check("closure_census", len(cands) == 59409 and len(active) == 2780
          and coeff_dist == {-1: 1357, 1: 1423} and covered_chains == 8344,
          "[COMPUTATION] exact closure (59409 candidates), sieve-active set (2780), "
          "coefficient distribution, covered chains (8344)")

    member_sets = [frozenset(m) for m in members]
    disagreements = 0
    for cuts in active:
        if explicit_subset_coefficient(frozenset(cuts), member_sets) != active[cuts][0]:
            disagreements += 1
    check("sieve_coefficient_agreement", disagreements == 0,
          "[COMPUTATION] hitting-mask sieve and explicit subset-union Mobius agree on every active term")
    stages["closure_and_sieve"] = time.process_time() - t0

    # Pair consistency census: which member pairs intersect (joint schedule
    # admissible).  This quantifies how much "no disjointness is assumed".
    t0 = time.process_time()
    ml = list(members)
    intersecting: list[list[int]] = []
    for i in range(len(ml)):
        for j in range(i + 1, len(ml)):
            if member_count(tuple(sorted(set(ml[i]) | set(ml[j])))) > 0:
                intersecting.append([i, j])
    check("member_pair_census", len(intersecting) == 2167,
          "[COMPUTATION] exactly 2167 of 81406 member pairs have an admissible joint schedule; "
          "the sieve, not pairwise disjointness, resolves all overlaps")
    stages["member_pair_census"] = time.process_time() - t0

    t0 = time.process_time()
    row_table: dict[tuple[tuple[int, int], ...], dict[str, int]] = {}
    for cuts in active:
        row_table[cuts] = d_row(cuts, Sched(cuts))
    stages["d_rows"] = time.process_time() - t0
    assert all(row_net(row_table[c]) != 0 for c in active), "some active term nets zero"

    # Old-41 regression: the full sieve projected onto the wave-13 sub-family
    # must reproduce the certified wave-13 union exactly.
    t0 = time.process_time()
    old_member_sets = {frozenset(sched2cuts(s)) for s in [s for _, s in old_specs]}
    old_h = sum(row_net(row_table[c]) for c in active if frozenset(c) in old_member_sets)
    base = (644233352324156721 + 333543290395947705 - 50032713330104219
            + 2 * 13710824278006626 + 2 * 2325336517026900)
    regression = base + old_h
    check("old41_regression", regression == WAVE13_UNION,
          f"[THEOREM] the wave-14 sieve projected onto the 41 wave-13 members "
          f"reproduces 507841195880595165555 exactly (delta {WAVE13_UNION - regression})")
    stages["regression"] = time.process_time() - t0

    overlap = {k: sum(c * row_table[cuts][k] for cuts, (c, _) in active.items())
              for k in ("H", "C", "L", "CL", "F2", "Q2", "F3", "Q3")}
    union = base + sum(c * row_net(row_table[cuts]) for cuts, (c, _) in active.items())
    naive_member = sum(row_net(row_table[m]) for m in members if m in active)
    naive_union = base + naive_member
    composite_sum = sum(c * row_net(row_table[cuts]) for cuts, (c, gnames) in active.items()
                       if len(gnames) >= 2)
    correction = naive_union - union
    pool = overlap["H"]
    member_pool = sum(member_count(m) for m in members)
    check("union_a35", union == 1972465461070186257835 and pool == 1972169603093083935789,
          "[THEOREM] exact finite sieve union gives a_35 >= 1972465461070186257835")
    check("overlap_breakdown", overlap == {"H": 1972169603093083935789, "C": 54142259738277507,
                                        "L": 629409190813341993, "CL": 50032678970365851,
                                        "F2": 13161781259087562, "Q2": 13161781259087562,
                                        "F3": 2057969889158220, "Q3": 2057969889158220},
          "[COMPUTATION] exact sieve aggregate of every intersection family")
    check("sieve_not_disjoint_assumption",
          naive_union == 7002669795485447916955 and correction == 5030204334415261659120
          and composite_sum == -correction and member_pool - pool == 5037180807639278388798,
          "[COMPUTATION] the disjoint-assumption value differs from the exact sieve union by the "
          "computed overlap correction; the sieve resolves 2167 intersecting pairs")

    m_value = C36 - union
    check("chain_M", m_value == 2939398390873631540302835,
          "[THEOREM] M = c_36 - a_35_lower from the first-backtrack chain")
    check("strict_improvement", union > WAVE13_UNION and m_value < C36 - WAVE13_UNION,
          "[THEOREM] strictly improves the wave-13 lower endpoint without using any K_c data")
    stages["aggregates"] = time.process_time() - t0

    lower, upper = atanh_root_interval(m_value)
    floor = floor40(lower)
    check("endpoint_floor", floor == "0.2122159753270231627267517174278577806020",
          "[COMPUTATION] 120-dps directed floor of atanh(M^(-1/36))")
    with localcontext() as context:
        context.prec = 220
        check("endpoint_rounding", Decimal(floor) <= Decimal(lower)
              and Decimal(upper) < Decimal(floor) + Decimal(1).scaleb(-REPORT_PLACES),
              "[COMPUTATION] floor <= interval lower <= true <= upper < floor + 10^-40")
    stages["endpoint"] = time.process_time() - t0

    active_rows = []
    for cuts in sorted(active):
        c, gnames = active[cuts]
        row = row_table[cuts]
        active_rows.append({
            "cuts": [[p, h] for p, h in cuts],
            "c": c,
            "members": sorted(gnames),
            "H": row["H"], "C": row["C"], "L": row["L"], "CL": row["CL"],
            "F2": row["F2"], "Q2": row["Q2"], "F3": row["F3"], "Q3": row["Q3"],
            "net": row_net(row),
        })
    member_rows = [{"name": name, "schedule": [c for pair in cuts for c in pair],
                   "count": member_count(cuts)} for cuts, name in members.items()]

    elapsed = time.process_time() - started
    resource_walls = [
        {
            "kind": "[COMPUTATION] observed",
            "description": "measured producer process time and peak resident set after all exact DPs",
            "seconds_process": round(elapsed, 6),
            "stages_process": {k: round(v, 6) for k, v in stages.items()},
            "peak_rss_bytes": peak_rss_bytes(),
        },
        {
            "kind": "[UNRESOLVED] preflight",
            "description": "declared input-size wall, not an observed lower/upper resource theorem",
            "profile_max_depth": 13,
            "members": 404,
            "closure_candidates": 59409,
            "active_terms": 2780,
            "member_pair_checks": 81406,
            "old41_separators": 820,
            "rss_budget_bytes": 8589934592,
        },
    ]
    payload = {
        "provenance": {
            "script": SCRIPT,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": ".venv/bin/python",
            "arithmetic": "[COMPUTATION] Python integers for all finite counts; mpmath.iv directed interval arithmetic at 120 dps only for the endpoint.",
            "mpmath_iv_dps": DPS,
            "elapsed_process_seconds": round(elapsed, 6),
            "stages_process": {k: round(v, 6) for k, v in stages.items()},
        },
        "data": {
            "classification": "[THEOREM] a_35 >= 1972465461070186257835 by an exact finite Mobius sieve union; [UNRESOLVED] this is not an all-size SAW-family theorem or a solution of the 3D Ising model.",
            "notation": {
                "c_r": "rooted r-step cubic-lattice SAW count",
                "a_r": "rooted r-step SAW count avoiding e_1=(1,0,0)",
                "h": "h(x,y,z)=-x+y+z; increments (-1,+1,+1,-1,+1,-1)",
                "M": "{-e_1,+e_2,+e_3}, indices (1,2,4)",
            },
            "external_source": {
                "status": "[EXTERNAL]",
                "path": "sources/fulltext/schram_barkema_bisseling2011.pdf",
                "sha256": source_hash,
                "citation": "Schram--Barkema--Bisseling (2011), J. Stat. Mech. P06019, Table I",
                "c_36": C36,
            },
            "positive_height_profile": {
                "status": "[COMPUTATION]",
                "method": "packed visited-set DFS of positive-height SAWs; strict-record terminal q is B(n,q)",
                "tail_T_0_to_13": tail,
                "blocks_B_n_q": {str(n): {str(q): count for q, count in sorted(blocks[n].items())}
                                for n in range(1, 14)},
                "parity_rule": "[LEMMA] every admissible rise q at block length l has q == l (mod 2)",
                "parity_verified": parity_ok,
            },
            "members": {
                "status": "[COMPUTATION]",
                "old41_count": 41,
                "two_block_admissible_count": 392,
                "total": len(members),
                "rows": member_rows,
            },
            "old41_disjointness": {
                "status": "[LEMMA]",
                "claim": "all 820 wave-13 schedule pairs have a recorded disjoint exact allowed-height slice",
                "pairwise_separator_rows": separator_rows,
            },
            "closure": {
                "status": "[COMPUTATION]",
                "method": "DP closure of admissible cut sequences over the member position universe",
                "candidates": len(cands),
                "g_size_census": {str(k): v for k, v in gsize_census.items()},
                "covered_chains": covered_chains,
                "active_terms": len(active),
                "composite_terms": sum(1 for _, (_, g) in active.items() if len(g) >= 2),
                "coefficient_distribution": coeff_dist,
                "sieve_method": "hitting-mask alternating sum; independently checked by the explicit subset-union Mobius sum (0 disagreements)",
            },
            "member_pair_overlaps": {
                "status": "[COMPUTATION]",
                "checked_pairs": 81406,
                "consistent_pairs": len(intersecting),
                "intersecting_pairs": intersecting,
                "statement": "[THEOREM] 2167 member pairs share an admissible joint schedule; no pairwise-disjointness assumption is made anywhere",
            },
            "intersection_families": {
                "status": "[THEOREM] inherited audited L, C, F2, Q2, F3, Q3 subfamilies of a_35",
                "L": 644233352324156721,
                "C": 333543290395947705,
                "C_cap_L": 50032713330104219,
                "F2_count": 13710824278006626,
                "Q2_count": 13710824278006626,
                "F3_count": 2325336517026900,
                "Q3_count": 2325336517026900,
                "base": base,
            },
            "active_terms": {
                "status": "[COMPUTATION]",
                "count": len(active_rows),
                "rows": active_rows,
            },
            "union_certificate": {
                "status": "[THEOREM]",
                "expression": "a_35 >= L + C - (C cap L) + F2 + Q2 + F3 + Q3 + sum_U c(U) (H_U - C cap H_U - L cap H_U + (C cap L) cap H_U - F2 cap H_U - Q2 cap H_U - F3 cap H_U - Q3 cap H_U)",
                "overlap_breakdown": overlap,
                "base": base,
                "naive_disjoint_union": naive_union,
                "overlap_correction": correction,
                "member_pool_sum": member_pool,
                "pool_union": pool,
                "a35_lower": union,
                "wave13_union": WAVE13_UNION,
                "delta_over_wave13": union - WAVE13_UNION,
                "scope": "[UNRESOLVED] finite r=35 certificate only; no all-r lower-family recursion is claimed.",
            },
            "chain": {
                "status": "[THEOREM]",
                "statement": "the previously proved first-backtrack inequality gives mu^n <= c_n - a_{n-1}; apply at n=36",
                "M": m_value,
                "mu_power_bound": "mu^36 <= 2939398390873631540302835",
            },
            "endpoint": {
                "status": "[THEOREM] via the SAW-domination proof in proofs/kc_bounds.md section A",
                "statement": "K_c >= atanh(M^(-1/36))",
                "interval_120dps": [lower, upper],
                "floor_40": floor,
                "rounding": "[COMPUTATION] {M} -> atanh(M^{-1/36}) -> directed floor: exact integer arithmetic down to the interval; the single rounding step floors the 120-dps interval lower bound at 40 decimal places with error term < 10^-40 (checked: floor <= lower <= true <= upper < floor + 10^-40)",
            },
            "resource_walls": resource_walls,
            "unresolved": [
                "[UNRESOLVED] the member family is the explicit 404-schedule subfamily (41 wave-13 plus 363 admissible two-block schedules), not an optimization over all cut schedules.",
                "[UNRESOLVED] no all-size strengthening of L_r or any Kesten-ratio theorem follows from this finite union.",
                "[UNRESOLVED] the finite certificate does not solve the three-dimensional Ising model.",
            ],
        },
        "checks": checks,
    }
    RESULT.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    print(f"\nwrote {RESULT}")
    print(f"elapsed process {elapsed:.1f}s; peak RSS {peak_rss_bytes()} bytes")
    print(f"a_35 lower bound: {union}")
    print(f"M: {m_value}")
    print(f"K_c lower floor (40): {floor}")


if __name__ == "__main__":
    try:
        main()
    except AssertionError as exc:
        print(f"CERTIFICATE FAILURE: {exc}")
        raise
