#!/usr/bin/env python3
"""Wave-11 finite path-family union certificate for the SAW K_c bound (e116).

Independently strengthen the wave-10 n=36 SAW certificate of
``proofs/saw_refine.md`` by taking the UNION of three exactly counted
sub-families of ``a_35`` (35-step self-avoiding walks avoiding ``e_1``)
with exact inclusion-exclusion, instead of the single height-slab family.

Exact content of this script:

1. [LEMMA] Height partition lemma.  With h(x,y,z)=-x+y+z every member of
   the slab family H has h = 7 exactly at step 11, while every defect-family
   member has at most one dh = -1 step, so its height at step 11 is 9 or 11.
   Hence H and the defect family L are DISJOINT - recomputed, not assumed
   (the prompt constants were treated as hypotheses and verified).
2. [COMPUTATION] The three finite families are recomputed from their
   definitions: |H| = B(11,7)^2 * T(13) with the two slab-block and
   half-space-tail enumerations; |L| = 3^35 + 4 G_35 + 2 D_35; and the
   coordinate-monotone family C = all 35-step words that are monotone in
   each coordinate with first step != +e_1, counted by a 27-state
   finite automaton, |C| = 333543290395947705.
3. [THEOREM] Exact overlaps.  C cap L = 3^35 + 4*(34*2^33) (each of the
   four corner-defect families contributes 34*2^33 coordinate-monotone
   words - their compatible alphabets have one letter per coordinate - and
   the fold families contribute none by sign reversal); C cap H is
   evaluated by an exact sign/height Moebius DP over the 27 orthant letter
   sets, |C cap H| = 157640944278888, validated against direct enumeration
   on three mini slab geometries.
4. [THEOREM] Union certificate.  a_35 >= |H| + |L| + |C| - |C cap L|
   - |C cap H| = 76401062626946721319, and the conservative sub-family
   E = C minus (M-words and the one-defect classes) reproduces the smaller
   fallback union 76401062592586982951.
5. [THEOREM] The first-backtrack submultiplicativity c_(m+n) <= c_m
   (c_n - a_(n-1)) is reproved and brute-force verified, iterated to
   mu^36 <= c_36 - a_35 <= M with M = 2941294455272074779839351, strictly
   below the wave-10 base c_36 - H_35.
6. [COMPUTATION] Directed mpmath interval arithmetic at 120 dps gives the
   new K_c lower endpoint floor
   0.2122120589214465859334619330429416874541 (40 places, rounded down).

No decision uses a floating-point critical-coupling benchmark.  Every
comparison is an exact integer comparison.
"""
from __future__ import annotations

from datetime import datetime, timezone
from decimal import Decimal, ROUND_FLOOR, localcontext
import hashlib
import json
import mpmath as mp
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

SCRIPT = "experiments/e116_saw_union.py"
RESULT_PATH = ROOT / "results" / "bounds" / "saw_union.json"
SAW_SOURCE = ROOT / "sources" / "fulltext" / "schram_barkema_bisseling2011.pdf"
DPS = 120
REPORT_PLACES = 40
C36 = 2941370856334701726560670
C36_SHA256 = "898d8333e8d70acd279e29f226ff108550b132463aa662f86df590b0ad72cb12"
INCUMBENT_K_LOWER = "0.2122119011661678393310862783954278184914"
WAVE10_FLOOR = "0.2122120570061123902353626109344172245745"
EXTERNAL_SAW_COUNTS = {
    30: 270569905525454674614,
    31: 1274191064726416905966,
    32: 5997359460809616886494,
    33: 28233744272563685150118,
    34: 132853629626823234210582,
    35: 625248129452557974777990,
    36: C36,
}
# In-repo stored exact series c_0..c_10 and a_0..a_8 (e05/kc_interval2
# enumerations, OEIS-anchored).  Independent anchors for this backtracking.
STORED_C_SMALL = (1, 6, 30, 150, 726, 3534, 16926, 81390, 387966, 1853886, 8809878)
STORED_A_SMALL = (1, 5, 25, 121, 589, 2821, 13565, 64661, 308981)

MAX_SERIES_DEPTH = 10
CM_BRUTE_DEPTH = 8          # brute force CM-vs-automaton agreement depth
L_BRUTE_DEPTH = 10          # defect-family brute-force verification depth
CHAIN_TOTAL_MAX = 10        # chain inequality brute-force verification depth
BLOCK_LEN, BLOCK_RISE, TAIL_LEN = 11, 7, 13
MINI_CONFIGS = ((5, 3, 4), (5, 3, 5), (5, 3, 6))

# Packed lattice points: 7 bits per signed coordinate, offset +32.
OFF = 32
STEPS = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))
DH = (-1, 1, 1, -1, 1, -1)            # h = -x+y+z
MONO_STEPS = (1, 2, 4)                # M = {-e1, +e2, +e3}
DEFECT_LETTERS = (0, 3, 5)            # +e1 (fold), -e2, -e3 (corner)


def pack(x: int, y: int, z: int) -> int:
    return ((x + OFF) << 14) | ((y + OFF) << 7) | (z + OFF)


def pack_delta(dx: int, dy: int, dz: int) -> int:
    # arithmetic sum: a bitwise OR is wrong when a lower field is negative
    return (dx << 14) + (dy << 7) + dz


PACKED_STEPS = tuple(pack_delta(*d) for d in STEPS)
FORBID_A = pack(1, 0, 0)

checks: list[dict[str, object]] = []


def _check(name: str, passed: bool, detail: str) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})
    if not passed:
        print(f"FAIL {name}: {detail}")
        raise AssertionError(name)
    print(f"PASS {name}: {detail}")


def _fingerprint(value: int) -> dict[str, object]:
    magnitude = abs(value)
    payload = magnitude.to_bytes(max(1, (magnitude.bit_length() + 7) // 8), "big")
    return {
        "sign": (value > 0) - (value < 0),
        "bit_length": magnitude.bit_length(),
        "sha256_magnitude_big_endian": hashlib.sha256(payload).hexdigest(),
    }


def _floor_decimal(value: str, places: int = REPORT_PLACES) -> str:
    with localcontext() as context:
        context.prec = places + 20
        quantum = Decimal(1).scaleb(-places)
        return format(Decimal(value).quantize(quantum, rounding=ROUND_FLOOR), "f")


def _interval_endpoints(interval) -> tuple[str, str]:
    text = str(interval).strip()
    lower, upper = text[1:-1].split(",", maxsplit=1)
    return lower.strip(), upper.strip()


def _atanh_neg_root_interval(magnitude: int, root: int) -> tuple[str, str]:
    """Directed-rounding enclosure of atanh(magnitude^(-1/root))."""
    previous = mp.iv.dps
    try:
        mp.iv.dps = DPS
        x = mp.iv.mpf(magnitude)
        v = mp.iv.exp(mp.iv.log(x) * (mp.iv.mpf(-1) / root))
        return _interval_endpoints(mp.iv.log((1 + v) / (1 - v)) / 2)
    finally:
        mp.iv.dps = previous


# ----------------------------------------------------------------------
# Enumerators (recursive, packed-int; the standalone verifier uses an
# iterative tuple-coordinate implementation).
# ----------------------------------------------------------------------

def enumerate_series(max_depth: int, forbidden: tuple[int, ...] = ()) -> list[int]:
    """Exact per-depth SAW counts avoiding `forbidden` vertices."""
    counts = [0] * (max_depth + 1)
    forbid = frozenset(forbidden)
    start = pack(0, 0, 0)

    def visit(point: int, depth: int) -> None:
        counts[depth] += 1
        if depth == max_depth:
            return
        for delta in PACKED_STEPS:
            nxt = point + delta
            if nxt in forbid or nxt in visited:
                continue
            visited.add(nxt)
            visit(nxt, depth + 1)
            visited.remove(nxt)

    visited = {start}
    visit(start, 0)
    return counts


def enumerate_step_words(max_depth: int) -> list[list[tuple[int, ...]]]:
    """SAW step-index words (avoiding e_1) up to max_depth, per length."""
    out: list[list[tuple[int, ...]]] = [[] for _ in range(max_depth + 1)]
    start = pack(0, 0, 0)
    word: list[int] = []

    def visit(point: int, depth: int) -> None:
        if depth >= 1:
            out[depth].append(tuple(word))
        if depth == max_depth:
            return
        for idx, delta in enumerate(PACKED_STEPS):
            nxt = point + delta
            if nxt == FORBID_A or nxt in visited:
                continue
            visited.add(nxt)
            word.append(idx)
            visit(nxt, depth + 1)
            word.pop()
            visited.remove(nxt)

    out[0].append(())  # the empty word is the unique depth-0 e_1-avoiding SAW
    visited = {start}
    visit(start, 0)
    return out


def count_height_block(length: int, rise: int) -> int:
    """SAWs of `length` steps from h=0 to h=rise, interior 0 < h < rise."""
    total = 0
    start = pack(0, 0, 0)

    def visit(point: int, height: int, depth: int) -> None:
        nonlocal total
        if depth == length:
            total += 1
            return
        for idx, delta in enumerate(PACKED_STEPS):
            nxt = point + delta
            if nxt in visited:
                continue
            nh = height + DH[idx]
            if depth + 1 == length:
                if nh != rise:
                    continue
            elif not 0 < nh < rise:
                continue
            visited.add(nxt)
            visit(nxt, nh, depth + 1)
            visited.remove(nxt)

    visited = {start}
    visit(start, 0, 0)
    return total


def count_halfspace_tail(length: int) -> int:
    """SAWs of `length` steps whose non-origin heights are positive."""
    total = 0
    start = pack(0, 0, 0)

    def visit(point: int, height: int, depth: int) -> None:
        nonlocal total
        if depth == length:
            total += 1
            return
        for idx, delta in enumerate(PACKED_STEPS):
            nxt = point + delta
            if nxt in visited:
                continue
            if height + DH[idx] <= 0:
                continue
            visited.add(nxt)
            visit(nxt, height + DH[idx], depth + 1)
            visited.remove(nxt)

    visited = {start}
    visit(start, 0, 0)
    return total


# ----------------------------------------------------------------------
# Defect family L: classifier and closed forms (independent of e92 code).
# ----------------------------------------------------------------------

def classify_l(word: tuple[int, ...]) -> str | None:
    """Family name of an e_1-avoiding SAW step word inside L, else None."""
    for s in word:
        if s not in MONO_STEPS:
            break
    else:
        return "mono"
    t = next(i for i, s in enumerate(word) if s not in MONO_STEPS)
    for s in word[t + 1:]:
        if s not in MONO_STEPS:
            return None
    bad = word[t]
    prev = word[t - 1] if t >= 1 else None
    nxt = word[t + 1] if t + 1 < len(word) else None
    if bad == 3:                                    # defect -e2
        if prev == 4 and nxt != 2:
            return "G1"
        if prev == 1 and nxt != 2:
            return "G3"
    elif bad == 5:                                  # defect -e3
        if prev == 2 and nxt != 4:
            return "G2"
        if prev == 1 and nxt != 4:
            return "G4"
    elif bad == 0:                                  # defect +e1 (fold)
        prev2 = word[t - 2] if t >= 2 else None
        if prev2 == 1 and prev == 2 and nxt != 1:
            return "D1"
        if prev2 == 1 and prev == 4 and nxt != 1:
            return "D2"
    return None


def defect_family_formulas(r: int) -> dict[str, int]:
    """Closed-form family counts L_r = 3^r + 4 G_r + 2 D_r."""

    def suffix_count(j: int) -> int:
        return 1 if j == 0 else 2 * 3 ** (j - 1)

    g = sum(3 ** k * suffix_count(r - 2 - k) for k in range(0, r - 1)) if r >= 2 else 0
    d = sum(3 ** k * suffix_count(r - 3 - k) for k in range(0, r - 2)) if r >= 3 else 0
    return {"mono": 3 ** r, "G_each": g, "D_each": d, "L_r": 3 ** r + 4 * g + 2 * d}


# ----------------------------------------------------------------------
# Coordinate-monotone family C: 27-state automaton.
# ----------------------------------------------------------------------

CM_AUTOMATON_DOC = (
    "State = (sx, sy, sz) with s in {unused, minus, plus} per coordinate (27 states). "
    "Transitions: append letter i (axis a = i//2, sign +1 if i even else -1) exactly when "
    "the coordinate slot is unused or matches; the slot is then fixed to that sign. "
    "The initial state is (unused, unused, unused) and the first letter is forbidden to be "
    "letter 0 (+e_1).  Accepting states: all, after exactly r letters.  A word is accepted "
    "iff it is coordinate-monotone and never visits e_1 (LEMMA: a coordinate-monotone word "
    "visits (1,0,0) iff its first letter is +e_1)."
)


def cm_automaton_count(r: int) -> int:
    from collections import defaultdict
    dp: dict[tuple[int, int, int], int] = defaultdict(int)
    dp[(0, 0, 0)] = 1
    for step in range(r):
        nd: dict[tuple[int, int, int], int] = defaultdict(int)
        for state, cnt in dp.items():
            for i in range(6):
                if i == 0 and step == 0:
                    continue
                ax = i // 2
                slot = 2 if i % 2 == 0 else 1
                cur = state[ax]
                if cur != 0 and cur != slot:
                    continue
                ns = list(state)
                ns[ax] = slot
                nd[tuple(ns)] += cnt
        dp = nd
    return sum(dp.values())


def is_cm_word(word: tuple[int, ...]) -> bool:
    slots = [0, 0, 0]
    for i in word:
        ax = i // 2
        slot = 2 if i % 2 == 0 else 1
        if slots[ax] != 0 and slots[ax] != slot:
            return False
        slots[ax] = slot
    return True


# ----------------------------------------------------------------------
# C cap L and the conservative removal classes, as finite automata.
# ----------------------------------------------------------------------

G_AUTOMATON_DOC = (
    "For each corner-defect family the coordinate-monotone members are exactly the words "
    "over the three-letter alphabet {-e_1, +e_other, defect} (one letter per coordinate, "
    "hence every such word is coordinate-monotone and avoids e_1) whose single defect "
    "letter is preceded by the prescribed letter.  Automaton state = (defect_seen, "
    "previous letter); counting words of length r gives (r-1)*2^(r-2).  The fold families "
    "contribute nothing: a fold word contains both -e_1 and +e_1, which reverses the "
    "x-sign slot, so no fold word is coordinate-monotone."
)


def g_family_cm_count(r: int, defect: int, prev_letter: int, free_letters: tuple[int, ...]) -> int:
    """Words of length r over free_letters + {defect} with exactly one defect,
    preceded by prev_letter, counted by the (defect_seen, prev) automaton."""
    if r < 2:
        return 0
    from collections import defaultdict
    dp: dict[tuple[int, int], int] = defaultdict(int)
    dp[(0, -1)] = 1
    for _ in range(r):
        nd: dict[tuple[int, int], int] = defaultdict(int)
        for (seen, prev), cnt in dp.items():
            for letter in free_letters:
                nd[(seen, letter)] += cnt
            if not seen and prev == prev_letter:
                nd[(1, defect)] += cnt
        dp = nd
    return sum(cnt for (seen, _prev), cnt in dp.items() if seen)


def one_defect_class_count(r: int, defect: int, free_letters: tuple[int, ...]) -> int:
    """Words of length r over free_letters + {defect} with exactly one defect
    (any position, no neighbourhood condition): the conservative class."""
    if r < 1:
        return 0
    from collections import defaultdict
    dp: dict[int, int] = defaultdict(int)
    dp[0] = 1
    for _ in range(r):
        nd: dict[int, int] = defaultdict(int)
        for seen, cnt in dp.items():
            for letter in free_letters:
                nd[seen] += cnt
            if not seen:
                nd[1] += cnt
        dp = nd
    return dp.get(1, 0)


# ----------------------------------------------------------------------
# C cap H: exact sign/height Moebius DP over orthant letter sets.
# ----------------------------------------------------------------------

MOBIUS_DOC = (
    "For an orthant letter set Q (at most one sign per coordinate; 27 sets incl. empty) "
    "let B1(Q) count 11-step words over Q from h=0 to h=7 with interior 0<h<7 and T(Q) "
    "count 13-step words over Q with all non-origin heights positive, both by an exact "
    "DP over integer heights (a word over Q is coordinate-monotone, hence self-avoiding, "
    "so no vertex bookkeeping is needed).  Because the slab family splits at steps 11 and "
    "22 and the second block is the height translate of the first, M(Q) = B1(Q)^2*T(Q) "
    "counts slab-family walks with steps in Q.  Counting each walk once by its used "
    "letter set gives |C cap H| = sum_Q (-1)^{e(Q)} M(Q) with e(Q) the number of "
    "coordinates empty in Q (Moebius inversion over the per-coordinate letter choices)."
)


def orthant_codes() -> list[tuple[int | None, int | None, int | None]]:
    from itertools import product
    return list(product((None, 0, 1), repeat=3))


def code_letters(code) -> tuple[int, ...]:
    out = []
    for ax, s in enumerate(code):
        if s is None:
            continue
        out.append(ax * 2 + (0 if s == 0 else 1))
    return tuple(out)


def block_dp(length: int, rise: int, letters: tuple[int, ...]) -> int:
    from collections import defaultdict
    dp: dict[int, int] = {0: 1}
    for step in range(length):
        nd: dict[int, int] = defaultdict(int)
        for h, cnt in dp.items():
            for i in letters:
                nh = h + DH[i]
                if step + 1 == length:
                    if nh != rise:
                        continue
                elif not 0 < nh < rise:
                    continue
                nd[nh] += cnt
        dp = nd
    return dp.get(rise, 0)


def tail_dp(length: int, letters: tuple[int, ...]) -> int:
    from collections import defaultdict
    dp: dict[int, int] = {0: 1}
    for _ in range(length):
        nd: dict[int, int] = defaultdict(int)
        for h, cnt in dp.items():
            for i in letters:
                nh = h + DH[i]
                if nh <= 0:
                    continue
                nd[nh] += cnt
        dp = nd
    return sum(dp.values())


def mobius_table(blk_len: int, blk_rise: int, tail_len: int) -> list[dict[str, object]]:
    rows = []
    for code in orthant_codes():
        letters = code_letters(code)
        if not letters:
            continue
        empty = sum(1 for s in code if s is None)
        bq = block_dp(blk_len, blk_rise, letters)
        tq = tail_dp(tail_len, letters)
        rows.append(
            {
                "letters": list(letters),
                "empty_coords": empty,
                "sign": (-1) ** empty,
                "B1_Q": bq,
                "T_Q": tq,
                "term": (-1) ** empty * bq * bq * tq,
            }
        )
    return rows


def c_cap_h_from_table(table: list[dict[str, object]]) -> int:
    return sum(int(row["term"]) for row in table)


def restricted_block_dfs(length: int, rise: int, letters: tuple[int, ...]) -> int:
    """Direct visited-set enumeration of blocks with steps restricted to Q."""
    allowed = frozenset(pack_delta(*STEPS[i]) for i in letters)
    total = 0
    start = pack(0, 0, 0)

    def visit(point: int, height: int, depth: int) -> None:
        nonlocal total
        if depth == length:
            total += 1
            return
        for idx, delta in enumerate(PACKED_STEPS):
            if delta not in allowed:
                continue
            nxt = point + delta
            if nxt in visited:
                continue
            nh = height + DH[idx]
            if depth + 1 == length:
                if nh != rise:
                    continue
            elif not 0 < nh < rise:
                continue
            visited.add(nxt)
            visit(nxt, nh, depth + 1)
            visited.remove(nxt)

    visited = {start}
    visit(start, 0, 0)
    return total


def restricted_tail_dfs(length: int, letters: tuple[int, ...]) -> int:
    allowed = frozenset(pack_delta(*STEPS[i]) for i in letters)
    total = 0
    start = pack(0, 0, 0)

    def visit(point: int, height: int, depth: int) -> None:
        nonlocal total
        if depth == length:
            total += 1
            return
        for idx, delta in enumerate(PACKED_STEPS):
            if delta not in allowed:
                continue
            nxt = point + delta
            if nxt in visited:
                continue
            if height + DH[idx] <= 0:
                continue
            visited.add(nxt)
            visit(nxt, height + DH[idx], depth + 1)
            visited.remove(nxt)

    visited = {start}
    visit(start, 0, 0)
    return total


# ----------------------------------------------------------------------
# Mini slab geometries: direct family enumeration (independent of the DP).
# ----------------------------------------------------------------------

def mini_family_count(blk_len: int, blk_rise: int, tail_len: int) -> int:
    """|H| for a mini slab geometry by direct DFS over the concatenated walks."""
    total = 0
    r = 2 * blk_len + tail_len
    start = pack(0, 0, 0)

    def visit(point: int, height: int, depth: int) -> None:
        nonlocal total
        if depth == r:
            total += 1
            return
        for idx, delta in enumerate(PACKED_STEPS):
            nxt = point + delta
            if nxt in visited:
                continue
            nh = height + DH[idx]
            if depth + 1 == blk_len:
                if nh != blk_rise:
                    continue
            elif depth < blk_len:
                if not 0 < nh < blk_rise:
                    continue
            elif depth + 1 == 2 * blk_len:
                if nh != 2 * blk_rise:
                    continue
            elif depth < 2 * blk_len:
                if not blk_rise < nh < 2 * blk_rise:
                    continue
            else:
                if nh <= 2 * blk_rise:
                    continue
            visited.add(nxt)
            visit(nxt, nh, depth + 1)
            visited.remove(nxt)

    visited = {start}
    visit(start, 0, 0)
    return total


def mini_cm_cap_h_direct(blk_len: int, blk_rise: int, tail_len: int) -> int:
    """Direct count of coordinate-monotone mini-slab walks (pruned DFS)."""
    total = 0
    r = 2 * blk_len + tail_len
    start = pack(0, 0, 0)
    slots = [0, 0, 0]

    def visit(point: int, height: int, depth: int, slots: list[int]) -> None:
        nonlocal total
        if depth == r:
            total += 1
            return
        for idx, delta in enumerate(PACKED_STEPS):
            nxt = point + delta
            if nxt in visited:
                continue
            nh = height + DH[idx]
            if depth + 1 == blk_len:
                if nh != blk_rise:
                    continue
            elif depth < blk_len:
                if not 0 < nh < blk_rise:
                    continue
            elif depth + 1 == 2 * blk_len:
                if nh != 2 * blk_rise:
                    continue
            elif depth < 2 * blk_len:
                if not blk_rise < nh < 2 * blk_rise:
                    continue
            else:
                if nh <= 2 * blk_rise:
                    continue
            ax = idx // 2
            slot = 2 if idx % 2 == 0 else 1
            if slots[ax] != 0 and slots[ax] != slot:
                continue
            saved = slots[ax]
            slots[ax] = slot
            visited.add(nxt)
            visit(nxt, nh, depth + 1, slots)
            visited.remove(nxt)
            slots[ax] = saved

    visited = {start}
    visit(start, 0, 0, slots)
    return total


def mini_h_cap_l_direct(blk_len: int, blk_rise: int, tail_len: int) -> int:
    """Direct count of mini-slab walks that are defect-family members.

    Pruned by the at-most-one-defect condition (any L member has at most
    one non-monotone letter), so the search is exhaustive for membership.
    """
    total = 0
    r = 2 * blk_len + tail_len
    start = pack(0, 0, 0)
    word: list[int] = []

    def visit(point: int, height: int, depth: int, defects: int) -> None:
        nonlocal total
        if depth == r:
            if classify_l(tuple(word)) is not None:
                total += 1
            return
        for idx, delta in enumerate(PACKED_STEPS):
            nxt = point + delta
            if nxt in visited:
                continue
            nh = height + DH[idx]
            if depth + 1 == blk_len:
                if nh != blk_rise:
                    continue
            elif depth < blk_len:
                if not 0 < nh < blk_rise:
                    continue
            elif depth + 1 == 2 * blk_len:
                if nh != 2 * blk_rise:
                    continue
            elif depth < 2 * blk_len:
                if not blk_rise < nh < 2 * blk_rise:
                    continue
            else:
                if nh <= 2 * blk_rise:
                    continue
            ndef = defects + (0 if idx in MONO_STEPS else 1)
            if ndef > 1:
                continue
            visited.add(nxt)
            word.append(idx)
            visit(nxt, nh, depth + 1, ndef)
            word.pop()
            visited.remove(nxt)

    visited = {start}
    visit(start, 0, 0, 0)
    return total


def main() -> None:
    started = time.monotonic()
    timings: dict[str, float] = {}

    # ------------------------------------------------------------------
    # External input audit and stored-series anchors
    # ------------------------------------------------------------------
    source_hash = hashlib.sha256(SAW_SOURCE.read_bytes()).hexdigest()
    _check(
        "external_saw_source_hash",
        source_hash == C36_SHA256,
        "cached Schram--Barkema--Bisseling PDF has the audited SHA-256",
    )

    # ------------------------------------------------------------------
    # Exact series (c and a) by fresh backtracking
    # ------------------------------------------------------------------
    t0 = time.monotonic()
    c_series = enumerate_series(MAX_SERIES_DEPTH)
    a_series = enumerate_series(MAX_SERIES_DEPTH, forbidden=(FORBID_A,))
    timings["series_enumeration"] = time.monotonic() - t0
    _check(
        "small_c_matches_stored_repo_series",
        tuple(c_series) == STORED_C_SMALL,
        "backtracking c_0..c_10 equals the OEIS-anchored stored series",
    )
    _check(
        "small_a_matches_stored_series",
        tuple(a_series[: len(STORED_A_SMALL)]) == STORED_A_SMALL,
        "backtracking a_0..a_8 equals the stored avoid-e_1 series",
    )

    # ------------------------------------------------------------------
    # Family H: slab blocks and half-space tail
    # ------------------------------------------------------------------
    t0 = time.monotonic()
    block_count = count_height_block(BLOCK_LEN, BLOCK_RISE)
    tail_count = count_halfspace_tail(TAIL_LEN)
    timings["height_enumeration"] = time.monotonic() - t0
    h35 = block_count * block_count * tail_count
    _check("slab_block_B_11_7", block_count == 729000,
           "direct finite enumeration of 11-step h-slab blocks ending at h=7")
    _check("halfspace_tail_T_13", tail_count == 142016661,
           "direct finite enumeration of 13-step positive-h half-space walks")
    _check("H35_product", h35 == 75473476338501000000,
           "H_35 = B(11,7)^2 * T(13) = 75473476338501000000")

    # ------------------------------------------------------------------
    # Family L: closed forms vs exhaustive classification
    # ------------------------------------------------------------------
    t0 = time.monotonic()
    words_by_depth = enumerate_step_words(L_BRUTE_DEPTH)
    l_rows: list[dict[str, object]] = []
    for r in range(L_BRUTE_DEPTH + 1):
        fam: dict[str, int] = {}
        cm_fam: dict[str, int] = {}
        for w in words_by_depth[r]:
            f = classify_l(w)
            if f is not None:
                fam[f] = fam.get(f, 0) + 1
                if is_cm_word(w):
                    cm_fam[f] = cm_fam.get(f, 0) + 1
        formulas = defect_family_formulas(r)
        exp_cm_g = (r - 1) * 2 ** (r - 2) if r >= 2 else 0
        ok = (
            fam.get("mono", 0) == formulas["mono"]
            and all(fam.get(f"G{i}", 0) == formulas["G_each"] for i in (1, 2, 3, 4))
            and all(fam.get(f"D{i}", 0) == formulas["D_each"] for i in (1, 2))
            and cm_fam.get("mono", 0) == 3 ** r
            and all(cm_fam.get(f"G{i}", 0) == exp_cm_g for i in (1, 2, 3, 4))
            and all(cm_fam.get(f"D{i}", 0) == 0 for i in (1, 2))
        )
        l_rows.append(
            {"r": r, "families": fam, "cm_families": cm_fam,
             "G_each_formula": formulas["G_each"], "cm_G_expected": exp_cm_g}
        )
        _check(f"defect_family_row_r{r}", ok,
               f"exhaustive classification at r={r}: family counts match formulas, "
               f"CM overlap is mono=3^r, each G=(r-1)2^(r-2), D=0")
    timings["defect_bruteforce"] = time.monotonic() - t0
    l35 = defect_family_formulas(35)["L_r"]
    _check("L35_formula", l35 == 644233352324156721,
           "L_35 = 3^35 + 4 G_35 + 2 D_35 = 644233352324156721")

    # Height-partition disjointness certificate (structural).
    _check(
        "height_partition_table",
        all(DH[i] == 1 for i in MONO_STEPS) and all(DH[i] == -1 for i in DEFECT_LETTERS),
        "every monotone letter has dh=+1 and every defect letter dh=-1, so an L-walk "
        "has height 11-2d in {9,11} at step 11 while H forces height exactly 7",
    )
    h_cap_l = 0

    # ------------------------------------------------------------------
    # Family C: automaton vs brute force
    # ------------------------------------------------------------------
    t0 = time.monotonic()
    cm_rows: list[dict[str, int]] = []
    for r in range(CM_BRUTE_DEPTH + 1):
        brute = sum(1 for w in words_by_depth[r] if is_cm_word(w))
        auto = cm_automaton_count(r)
        cm_rows.append({"r": r, "bruteforce_cm_avoiding_e1": brute, "automaton": auto})
        _check(f"cm_automaton_r{r}", brute == auto,
               f"CM words avoiding e_1 by exhaustive SAW enumeration = automaton count at r={r}")
    c35 = cm_automaton_count(35)
    _check("C35_automaton", c35 == 333543290395947705,
           "|C| = 333543290395947705 by the 27-state automaton")
    timings["cm_checks"] = time.monotonic() - t0

    # ------------------------------------------------------------------
    # C cap L by per-family automata
    # ------------------------------------------------------------------
    mono_cm = 3 ** 35
    g_specs = (
        ("G1", 3, 4, (1, 4)),   # defect -e2 preceded by +e3, free letters {-e1,+e3}
        ("G2", 5, 2, (1, 2)),   # defect -e3 preceded by +e2, free letters {-e1,+e2}
        ("G3", 3, 1, (1, 4)),   # defect -e2 preceded by -e1
        ("G4", 5, 1, (1, 2)),   # defect -e3 preceded by -e1
    )
    g_cm_counts = {}
    for name, defect, prev_letter, free in g_specs:
        dp_val = g_family_cm_count(35, defect, prev_letter, free)
        closed = 34 * 2 ** 33
        g_cm_counts[name] = dp_val
        _check(f"cm_cap_{name}_dp_matches_closed_form", dp_val == closed,
               f"{name} coordinate-monotone members: automaton = (r-1)2^(r-2) = 34*2^33 = {dp_val}")
    d_cm = 0
    _check("cm_cap_D_zero_structural", True,
           "fold words contain both -e_1 and +e_1 (sign reversal), so no fold word is "
           "coordinate-monotone; brute force confirmed D=0 for r<=10 above")
    c_cap_l = mono_cm + sum(g_cm_counts.values()) + d_cm
    _check("C_cap_L_value", c_cap_l == 50032713330104219,
           "|C cap L| = 3^35 + 4*34*2^33 = 50032713330104219")

    # ------------------------------------------------------------------
    # C cap H: Moebius DP + independent verifications
    # ------------------------------------------------------------------
    t0 = time.monotonic()
    table = mobius_table(BLOCK_LEN, BLOCK_RISE, TAIL_LEN)
    c_cap_h = c_cap_h_from_table(table)
    _check("C_cap_H_mobius", c_cap_h == 157640944278888,
           "|C cap H| = 157640944278888 by the orthant Moebius DP")

    # per-Q DP vs direct restricted visited-set enumeration at full scale
    for row in table:
        letters = tuple(row["letters"])
        bdfs = restricted_block_dfs(BLOCK_LEN, BLOCK_RISE, letters)
        tdfs = restricted_tail_dfs(TAIL_LEN, letters)
        if bdfs != int(row["B1_Q"]) or tdfs != int(row["T_Q"]):
            _check(f"restricted_dfs_Q_{''.join(map(str, letters))}", False,
                   "restricted DFS disagrees with the height DP")
            break
    else:
        _check("restricted_dfs_all_orthants", True,
               "for all 26 nonempty orthants the restricted visited-set enumerations "
               "reproduce B1(Q) and T(Q) at the full (11,7,13) scale")

    mini_rows = []
    for (bl, br, tl) in MINI_CONFIGS:
        fam_n = mini_family_count(bl, br, tl)
        prod = None
        mb = c_cap_h_from_table(mobius_table(bl, br, tl))
        direct = mini_cm_cap_h_direct(bl, br, tl)
        hl = mini_h_cap_l_direct(bl, br, tl)
        mini_rows.append(
            {"config": [bl, br, tl], "|H_mini|": fam_n, "mobius_C_cap_H": mb,
             "direct_C_cap_H": direct, "direct_H_cap_L": hl}
        )
        _check(f"mini_mobius_{bl}_{br}_{tl}", mb == direct,
               f"mini slab ({bl},{br},{tl}): direct CM count {direct} = Moebius DP {mb}")
        _check(f"mini_H_cap_L_{bl}_{br}_{tl}", hl == 0,
               f"mini slab ({bl},{br},{tl}): no member is a defect-family walk")
    _check("mini_family_count_consistency",
           mini_rows[0]["|H_mini|"] == 108 * 108 * 171
           and mini_rows[1]["|H_mini|"] == 108 * 108 * 837
           and mini_rows[2]["|H_mini|"] == 108 * 108 * 3411,
           "mini family sizes equal B(5,3)^2*T(tl) with B(5,3)=108, T(4,5,6)=171,837,3411")
    timings["mobius_checks"] = time.monotonic() - t0

    # ------------------------------------------------------------------
    # Union certificate (exact inclusion-exclusion) and conservative E
    # ------------------------------------------------------------------
    h_union_l = h35 + l35
    _check("H_union_L_value", h_union_l == 76117709690825156721,
           "H cap L = empty, so |H union L| = H_35 + L_35 = 76117709690825156721")
    union = h_union_l + c35 - c_cap_l - c_cap_h
    _check("union_value", union == 76401062626946721319,
           "|H cup L cup C| = 76401062626946721319 by three-set inclusion-exclusion "
           "(all pairwise overlaps exact, triple overlap inside H cap L = empty)")

    r1 = one_defect_class_count(35, 3, (1, 4))
    r2 = one_defect_class_count(35, 5, (1, 2))
    _check("removal_class_sizes", r1 == 35 * 2 ** 34 and r2 == 35 * 2 ** 34,
           "each one-defect removal class has 35*2^34 words")
    e35 = c35 - 3 ** 35 - r1 - r2
    _check("E35_value", e35 == 283510542706105118,
           "|E| = |C| - 3^35 - 2*35*2^34 = 283510542706105118")
    e_net = e35 - c_cap_h
    union_e = h_union_l + e_net
    _check("E_route_union", e_net == 283352901761826230
           and union_e == 76401062592586982951,
           "conservative E route: net 283352901761826230, union 76401062592586982951 "
           "(the supplied prompt constants, now derived)")
    _check("E_strictly_weaker", union_e < union and union - union_e == 4 * 2 ** 33,
           "the E route is weaker than the exact C route by exactly 4*2^33")

    # brute-force the E cap L = empty claim at small r
    e_ok = True
    for r in range(L_BRUTE_DEPTH + 1):
        rc1 = one_defect_class_count(r, 3, (1, 4))
        rc2 = one_defect_class_count(r, 5, (1, 2))
        ec = cm_automaton_count(r) - 3 ** r - rc1 - rc2
        for w in words_by_depth[r]:
            if is_cm_word(w) and classify_l(w) is not None:
                # every CM L-member must lie in mono or a removal class
                mono_word = all(s in MONO_STEPS for s in w)
                in_r1 = w.count(3) == 1 and all(s in (1, 4) for s in w if s != 3)
                in_r2 = w.count(5) == 1 and all(s in (1, 2) for s in w if s != 5)
                if not (mono_word or in_r1 or in_r2):
                    e_ok = False
        if ec < 0:
            e_ok = False
    _check("E_cap_L_empty_bruteforce", e_ok,
           "for r<=10 every coordinate-monotone defect-family word lies in the mono or "
           "one-defect removal classes, certifying E cap L = empty at those depths")

    # ------------------------------------------------------------------
    # Chain theorem: brute-force verification and iteration metadata
    # ------------------------------------------------------------------
    chain_rows = []
    for total in range(2, CHAIN_TOTAL_MAX + 1):
        for n in range(1, total):
            m = total - n
            lhs = c_series[total]
            rhs = c_series[m] * (c_series[n] - a_series[n - 1])
            chain_rows.append({"m": m, "n": n, "lhs": lhs, "rhs": rhs, "passed": lhs <= rhs})
            _check(f"first_backtrack_m{m}_n{n}", lhs <= rhs,
                   "c_(m+n) <= c_m (c_n - a_(n-1))")
    small_chain = []
    for n in (2, 3):
        x_n = c_series[n] - a_series[n - 1]
        for k in range(2, MAX_SERIES_DEPTH // n + 1):
            if k * n > MAX_SERIES_DEPTH:
                continue
            lhs = c_series[k * n]
            rhs = c_series[n] * x_n ** (k - 1)
            small_chain.append({"n": n, "k": k, "lhs": lhs, "rhs": rhs, "passed": lhs <= rhs})
            _check(f"chain_n{n}_k{k}", lhs <= rhs, "c_(kn) <= c_n (c_n - a_(n-1))^(k-1)")

    m_value = C36 - union
    m_e92 = C36 - h35
    m_e = C36 - union_e
    _check("M_value", m_value == 2941294455272074779839351,
           "M = c_36 - |H cup L cup C| = 2941294455272074779839351")
    _check("M_strictly_below_wave10", m_value < m_e92,
           f"M = {m_value} < c_36 - H_35 = {m_e92}: the union strictly strengthens the wave-10 base")
    _check("M_E_between", m_e92 > m_e > m_value,
           "the conservative E-route base sits strictly between the wave-10 base and M")
    _check("chain_instance_beats_incumbent", m_value ** 36 < C36 ** 36,
           "mu^36 <= M < c_36 implies strict improvement over mu <= c_36^(1/36)")
    improvement_int = m_value ** 36 - C36 ** 36

    # honest finite selection record (integer comparisons only)
    external_lower = {n: defect_family_formulas(n - 1)["L_r"] for n in EXTERNAL_SAW_COUNTS if n != 36}
    external_lower[36] = union
    selection_rows = []
    for n in sorted(EXTERNAL_SAW_COUNTS):
        base = EXTERNAL_SAW_COUNTS[n] - external_lower[n]
        beats = base ** 36 < m_value ** n
        selection_rows.append(
            {"n": n, "c_n": EXTERNAL_SAW_COUNTS[n],
             "lower_bound_on_a_(n-1)": external_lower[n],
             "name": "union_35" if n == 36 else f"L_{n-1}",
             "base": base, "base_beats_our_certificate": beats}
        )
        _check(f"selection_n{n}_honest", n == 36 or not beats,
               f"at n={n} the recorded finite lower family does not beat the n=36 union "
               f"certificate (exact integer comparison)")

    # ------------------------------------------------------------------
    # High-precision endpoint
    # ------------------------------------------------------------------
    new_lo, new_hi = _atanh_neg_root_interval(m_value, 36)
    old_lo, old_hi = _atanh_neg_root_interval(m_e92, 36)
    inc_lo, inc_hi = _atanh_neg_root_interval(C36, 36)
    e_lo, e_hi = _atanh_neg_root_interval(m_e, 36)
    new_floor = _floor_decimal(new_lo)
    old_floor = _floor_decimal(old_lo)
    inc_floor = _floor_decimal(inc_lo)
    e_floor = _floor_decimal(e_lo)
    _check("incumbent_endpoint_reproduced", inc_floor == INCUMBENT_K_LOWER,
           "atanh(c_36^(-1/36)) reproduces the pre-wave-10 incumbent floor")
    _check("wave10_endpoint_reproduced", old_floor == WAVE10_FLOOR,
           "atanh((c_36-H_35)^(-1/36)) reproduces the wave-10 floor")
    _check("endpoint_floor_value", new_floor == "0.2122120589214465859334619330429416874541",
           "new 40-place downward floor of atanh(M^(-1/36))")
    with localcontext() as context:
        context.prec = REPORT_PLACES + 80
        _check(
            "directed_floor_is_true_lower_bound",
            Decimal(new_floor) <= Decimal(new_lo)
            and Decimal(new_hi) < Decimal(new_floor) + Decimal(1).scaleb(-REPORT_PLACES),
            "floor(x) <= lower endpoint < upper endpoint < floor + 10^-40 at 120 dps",
        )
        _check(
            "strict_improvement_chain",
            Decimal(inc_floor) < Decimal(old_floor) < Decimal(e_floor) < Decimal(new_floor),
            "incumbent < wave-10 < E-route < union endpoint (each strict)",
        )

    elapsed = time.monotonic() - started

    payload = {
        "provenance": {
            "script": SCRIPT,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": ".venv/bin/python",
            "arithmetic": (
                "All decisions are Python integer comparisons; mpmath.iv with directed "
                "rounding at 120 dps formats endpoints only."
            ),
            "mpmath_iv_dps": DPS,
            "elapsed_seconds": elapsed,
            "timings": timings,
        },
        "data": {
            "classification": (
                "[THEOREM] the three-family exact union a_35 >= |H cup L cup C| "
                "= 76401062626946721319 with exactly computed overlaps strictly strengthens "
                "the wave-10 height-slab certificate, giving mu^36 <= M and a new K_c lower "
                "endpoint; [LEMMA] H cap L = empty by the step-11 height partition; the "
                "conservative E-route union reproducing the supplied prompt constants is "
                "derived as a strict fallback; the union bound is a finite n=35 family "
                "statement, NOT an all-size family theorem; mu^2 <= c_32/c_30 remains open."
            ),
            "notation": {
                "c_r": "rooted r-step self-avoiding walks on Z^3",
                "a_r": "rooted r-step self-avoiding walks avoiding the vertex e_1",
                "h": "h(x,y,z) = -x+y+z with per-letter dh = (-1,1,1,-1,1,-1)",
                "letters": "+e1,-e1,+e2,-e2,+e3,-e3 = indices 0..5",
                "M_set": "monotone letter set {-e1,+e2,+e3} = indices (1,2,4), all dh=+1",
                "families": {
                    "H": "concatenations block*block*tail: 11-step SAWs h:0->7 interior "
                         "(0,7), translated second block h:7->14 interior (7,14), 13-step "
                         "tail with all post-junction heights > 14",
                    "L": "mono plus four one-defect corner families plus two fold families "
                         "(wave-10 defect family, definitions re-derived here)",
                    "C": "35-step words monotone in each coordinate with first letter "
                         "!= +e_1 (equivalently: coordinate-monotone 35-SAWs avoiding e_1)",
                    "E": "C minus the M-words minus the class {exactly one -e2, others in "
                         "{-e1,+e3}} minus the z-analogue (conservative sub-family)",
                },
            },
            "external_source": {
                "path": "sources/fulltext/schram_barkema_bisseling2011.pdf",
                "sha256": source_hash,
                "citation": "Schram, Barkema, Bisseling (2011), J. Stat. Mech. P06019, Table I",
                "c_36": C36,
                "external_c_30_to_36": {str(k): v for k, v in EXTERNAL_SAW_COUNTS.items()},
            },
            "series": {"c_0_to_10": c_series, "a_0_to_10": a_series},
            "family_H": {
                "status": "[COMPUTATION]",
                "block_len": BLOCK_LEN, "block_rise": BLOCK_RISE, "tail_len": TAIL_LEN,
                "B_11_7": block_count, "T_13": tail_count, "H_35": h35,
                "injectivity": (
                    "blocks/tail have disjoint height bands and meet only at the splice "
                    "vertices at steps 11 and 22, so the concatenation is self-avoiding and "
                    "the triple is recovered from the walk by its fixed split points; "
                    "h(e_1)=-1<1<=every non-origin height, so H avoids e_1"
                ),
            },
            "family_L": {
                "status": "[THEOREM] a_r >= L_r (wave-10, reproved here)",
                "formulas": "s(0)=1, s(j)=2*3^(j-1); G_r=sum_k 3^k s(r-2-k); "
                            "D_r=sum_k 3^k s(r-3-k); L_r=3^r+4G_r+2D_r",
                "L_35": l35,
                "L_29": defect_family_formulas(29)["L_r"],
                "bruteforce_rows": l_rows,
            },
            "family_C": {
                "status": "[COMPUTATION]",
                "automaton": CM_AUTOMATON_DOC,
                "C_35": c35,
                "lemma_e1": (
                    "a coordinate-monotone word visits e_1=(1,0,0) iff its first letter is "
                    "+e_1: reaching x=1 needs the x-sign + and exactly one +e_1 step, and "
                    "y=z=0 needs no y/z step before it, so the prefix is the single first "
                    "step; brute-force confirmed for r<=8"
                ),
                "bruteforce_rows": cm_rows,
            },
            "overlaps": {
                "H_cap_L": {
                    "status": "[LEMMA]",
                    "value": h_cap_l,
                    "proof": (
                        "every H-walk has height exactly 7 at step 11; every L-walk has at "
                        "most one dh=-1 letter (mono letters all have dh=+1, each family has "
                        "exactly one defect letter with dh=-1), so its step-11 height is "
                        "11-2d with d in {0,1} = {11,9}; 7 is in neither"
                    ),
                    "mini_verification": mini_rows,
                },
                "C_cap_L": {
                    "status": "[THEOREM]",
                    "value": c_cap_l,
                    "mono_cap_C": mono_cm,
                    "G_each_cap_C": {k: v for k, v in g_cm_counts.items()},
                    "D_cap_C": d_cm,
                    "automata": G_AUTOMATON_DOC,
                    "closed_form": "3^35 + 4*(34*2^33)",
                    "bruteforce_rows": [
                        {"r": row["r"], "cm_families": row["cm_families"],
                         "cm_G_expected": row["cm_G_expected"]} for row in l_rows
                    ],
                },
                "C_cap_H": {
                    "status": "[COMPUTATION]",
                    "value": c_cap_h,
                    "method": MOBIUS_DOC,
                    "mobius_table": table,
                    "restricted_dfs_crosscheck": (
                        "all 26 nonempty orthants re-enumerated with visited-set DFS at "
                        "the full (11,7,13) scale agree with the height DP"
                    ),
                    "mini_verification": mini_rows,
                },
                "triple": {"value": 0, "reason": "H cap L = empty contains the triple overlap"},
            },
            "union_certificate": {
                "status": "[THEOREM]",
                "expression": "|H| + |L| + |C| - |C cap L| - |C cap H| - |H cap L|",
                "value": union,
                "components": {
                    "H_35": h35, "L_35": l35, "H_union_L": h_union_l, "C_35": c35,
                    "C_cap_L": c_cap_l, "C_cap_H": c_cap_h, "H_cap_L": h_cap_l,
                },
                "scope": "a finite n=35 sub-family statement, not an all-size family",
            },
            "conservative_E_route": {
                "status": "[THEOREM] (strictly weaker fallback)",
                "R1_size": r1, "R2_size": r2, "E_35": e35,
                "E_cap_H": c_cap_h,
                "E_cap_H_reason": (
                    "M-words have no dh=-1 letter and each removal-class word exactly one, "
                    "so their step-11 heights lie in {11,9}, never 7; hence "
                    "E cap H = C cap H"
                ),
                "E_cap_L": 0,
                "E_cap_L_reason": (
                    "every coordinate-monotone L-word lies in the mono family or one of the "
                    "two removal classes (fold words are never coordinate-monotone), so "
                    "removing those classes empties E cap L; brute-force r<=10"
                ),
                "H_union_L": h_union_l,
                "E_net": e_net,
                "union_E": union_e,
                "gap_to_exact": union - union_e,
            },
            "chain_theorem": {
                "status": "[THEOREM]",
                "statement": (
                    "for all m,n >= 1: c_(m+n) <= c_m (c_n - a_(n-1)); iterating at fixed n "
                    "gives c_(kn) <= c_n (c_n - a_(n-1))^(k-1) and, by Fekete's lemma along "
                    "multiples of n, mu^n <= c_n - a_(n-1)"
                ),
                "proof_sketch": (
                    "split an (m+n)-SAW after step m into an m-SAW prefix and a translated "
                    "n-SAW suffix; valid pairs biject with (m+n)-SAWs.  For fixed prefix "
                    "omega with last step sigma, the n-SAWs starting with the reversal -sigma "
                    "followed by any (n-1)-SAW from omega_(m-1) avoiding omega_m form an "
                    "injected invalid family of size a_(n-1) (cubic symmetry fixes the "
                    "avoided neighbour).  Hence at most c_n - a_(n-1) suffixes survive per "
                    "prefix."
                ),
                "bruteforce_rows": chain_rows,
                "small_chain_rows": small_chain,
                "n36_instance": {
                    "a_35_lower": union,
                    "M": m_value,
                    "M_expression": "c_36 - |H cup L cup C|",
                    "M_wave10": m_e92,
                    "M_E_route": m_e,
                    "improvement_integer": {
                        "expression": "M^36 - c_36^36",
                        "fingerprint": _fingerprint(improvement_int),
                    },
                },
            },
            "endpoint": {
                "status": "[THEOREM] via kc_bounds SAW domination",
                "statement": "K_c >= atanh(M^(-1/36))",
                "new_interval_120dps": [new_lo, new_hi],
                "wave10_interval_120dps": [old_lo, old_hi],
                "incumbent_interval_120dps": [inc_lo, inc_hi],
                "E_route_interval_120dps": [e_lo, e_hi],
                "new_floor_40": new_floor,
                "wave10_floor_40": old_floor,
                "E_route_floor_40": e_floor,
                "incumbent_floor_40": inc_floor,
                "strictly_improves": True,
            },
            "selection_record": {
                "rows": selection_rows,
                "discipline": (
                    "n=36 is the largest exactly counted external c_n and the only one with "
                    "a union lower family; all comparisons are integer comparisons; no "
                    "floating critical-coupling benchmark was consulted anywhere"
                ),
            },
            "unresolved": [
                {
                    "status": "[UNRESOLVED]",
                    "statement": "the union family exists only at n=35; no all-r lower "
                                 "family strengthening L_r is constructed here",
                },
                {
                    "status": "[UNRESOLVED]",
                    "statement": "mu^2 <= c_32/c_30 would still need "
                                 "a_29 >= c_30 - floor(c_32/c_30) = 270569905525454674592 "
                                 "against L_29 = 741377533262625",
                },
            ],
        },
        "checks": checks,
    }

    RESULT_PATH.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    print(f"\nwrote {RESULT_PATH}")
    print(f"elapsed {elapsed:.1f}s")
    print(f"union a_35 lower bound: {union}")
    print(f"M = {m_value}")
    print(f"new lower endpoint (floor 40): {new_floor}")
    print(f"wave-10 endpoint (floor 40):   {old_floor}")


if __name__ == "__main__":
    try:
        main()
    except AssertionError as exc:
        print(f"CERTIFICATE FAILURE: {exc}")
        raise
