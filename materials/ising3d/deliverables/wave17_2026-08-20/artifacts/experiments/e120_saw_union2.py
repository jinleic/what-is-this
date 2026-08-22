#!/usr/bin/env python3
"""Exact wave-12 SAW-family union certificate at n=36.

[THEOREM] This program constructs a finite union of e_1-avoiding 35-step
self-avoiding-walk families.  Its new part is a disjoint 5-by-5 grid of
height slabs, plus two prescribed two-defect monotone languages.  Every
finite overlap used in inclusion--exclusion is counted by an exact finite
DP/Möbius calculation; no critical-coupling benchmark is consulted.

The standalone verifier is intentionally clean-room and recomputes the
positive-height profile by tuple-coordinate iterative DFS at 150 dps for the
endpoint enclosure.  Run this producer before the verifier:

    PYTHONPATH=src .venv/bin/python experiments/e120_saw_union2.py
    PYTHONPATH=src .venv/bin/python tests/test_saw_union2.py
"""
from __future__ import annotations

from collections import defaultdict
from datetime import datetime, timezone
from decimal import Decimal, ROUND_FLOOR, localcontext
from itertools import product
import hashlib
import json
import mpmath as mp
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

SCRIPT = "experiments/e120_saw_union2.py"
RESULT = ROOT / "results" / "bounds" / "saw_union2.json"
SAW_SOURCE = ROOT / "sources" / "fulltext" / "schram_barkema_bisseling2011.pdf"

C36 = 2941370856334701726560670
C36_SHA256 = "898d8333e8d70acd279e29f226ff108550b132463aa662f86df590b0ad72cb12"
DPS = 120
REPORT_PLACES = 40
R = 35
M_LETTERS = (1, 2, 4)  # {-e_1,+e_2,+e_3}; each has Δh=+1.
STEPS = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))
DH = (-1, 1, 1, -1, 1, -1)
OFF = 32

checks: list[dict[str, object]] = []


def check(name: str, passed: bool, detail: str) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})
    if not passed:
        print(f"FAIL {name}: {detail}")
        raise AssertionError(name)
    print(f"PASS {name}: {detail}")


def pack(x: int, y: int, z: int) -> int:
    return ((x + OFF) << 14) | ((y + OFF) << 7) | (z + OFF)


def pack_delta(dx: int, dy: int, dz: int) -> int:
    return (dx << 14) + (dy << 7) + dz


PACKED_STEPS = tuple(pack_delta(*step) for step in STEPS)


# ---------------------------------------------------------------------------
# Exact positive-height profile.
# ---------------------------------------------------------------------------


def positive_height_profile(max_depth: int) -> tuple[list[int], list[dict[int, int]]]:
    """Return T(n) and B(n,q) for all n <= ``max_depth`` by visited-set DFS.

    [LEMMA] A positive-height SAW ending at q is a B(n,q) block iff q is a
    strict record height along the path: every earlier non-origin height is
    then in (0,q), and conversely.  Thus one exact DFS produces all profile
    entries needed by the declared bounded split search.
    """
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
# Height schedules and exact coordinate-monotone intersections.
# ---------------------------------------------------------------------------


def schedule_allows(height: int, position: int, schedule: tuple[int, int, int, int, int]) -> bool:
    """Height rule for two blocks (l1,q1),(l2,q2), then a positive tail."""
    l1, q1, l2, q2, _tail = schedule
    if position < l1:
        return 0 < height < q1
    if position == l1:
        return height == q1
    if position < l1 + l2:
        return q1 < height < q1 + q2
    if position == l1 + l2:
        return height == q1 + q2
    return height > q1 + q2


def block_dp(length: int, rise: int, letters: tuple[int, ...]) -> int:
    """Exact word DP for a coordinate-monotone block over a signed orthant."""
    dp: dict[int, int] = {0: 1}
    for position in range(1, length + 1):
        next_dp: dict[int, int] = defaultdict(int)
        for height, count in dp.items():
            for letter in letters:
                next_height = height + DH[letter]
                if position == length:
                    if next_height != rise:
                        continue
                elif not 0 < next_height < rise:
                    continue
                next_dp[next_height] += count
        dp = next_dp
    return dp.get(rise, 0)


def tail_dp(length: int, letters: tuple[int, ...]) -> int:
    """Exact word DP for a coordinate-monotone positive-height tail."""
    dp: dict[int, int] = {0: 1}
    for _ in range(length):
        next_dp: dict[int, int] = defaultdict(int)
        for height, count in dp.items():
            for letter in letters:
                next_height = height + DH[letter]
                if next_height > 0:
                    next_dp[next_height] += count
        dp = next_dp
    return sum(dp.values())


def orthant_codes() -> list[tuple[int | None, int | None, int | None]]:
    return list(product((None, 0, 1), repeat=3))


def code_letters(code: tuple[int | None, int | None, int | None]) -> tuple[int, ...]:
    return tuple(axis * 2 + (0 if sign == 0 else 1)
                 for axis, sign in enumerate(code) if sign is not None)


def cm_automaton_count(length: int) -> int:
    """Count C with the finite 27-state signed-coordinate automaton."""
    dp: dict[tuple[int, int, int], int] = {(0, 0, 0): 1}
    for position in range(length):
        next_dp: dict[tuple[int, int, int], int] = defaultdict(int)
        for slots, count in dp.items():
            for letter in range(6):
                if position == 0 and letter == 0:
                    continue
                axis = letter // 2
                sign = 2 if letter % 2 == 0 else 1
                if slots[axis] not in (0, sign):
                    continue
                next_slots = list(slots)
                next_slots[axis] = sign
                next_dp[tuple(next_slots)] += count
        dp = next_dp
    return sum(dp.values())


def cm_overlap_mobius(schedule: tuple[int, int, int, int, int]) -> tuple[int, list[dict[str, object]]]:
    """[THEOREM] C∩H by signed-orthant Möbius inversion."""
    l1, q1, l2, q2, tail_length = schedule
    rows: list[dict[str, object]] = []
    for code in orthant_codes():
        letters = code_letters(code)
        if not letters:
            continue
        empty = sum(sign is None for sign in code)
        b1 = block_dp(l1, q1, letters)
        b2 = block_dp(l2, q2, letters)
        tail = tail_dp(tail_length, letters)
        term = (-1) ** empty * b1 * b2 * tail
        rows.append({
            "code": [sign for sign in code],
            "letters": list(letters),
            "empty_coordinates": empty,
            "sign": (-1) ** empty,
            "B1": b1,
            "B2": b2,
            "T": tail,
            "term": term,
        })
    return sum(int(row["term"]) for row in rows), rows


# ---------------------------------------------------------------------------
# Exact L-language intersections with a height schedule.
# ---------------------------------------------------------------------------


def language_schedule_count(letters_by_position: list[tuple[int, ...]],
                            schedule: tuple[int, int, int, int, int],
                            coordinate_monotone: bool = False) -> int:
    """Exact finite language DP, optionally retaining three sign slots."""
    if coordinate_monotone:
        dp: dict[tuple[int, tuple[int, int, int]], int] = {(0, (0, 0, 0)): 1}
        for position, letters in enumerate(letters_by_position, start=1):
            next_dp: dict[tuple[int, tuple[int, int, int]], int] = defaultdict(int)
            for (height, slots), count in dp.items():
                for letter in letters:
                    next_height = height + DH[letter]
                    if not schedule_allows(next_height, position, schedule):
                        continue
                    axis = letter // 2
                    sign = 2 if letter % 2 == 0 else 1
                    if slots[axis] not in (0, sign):
                        continue
                    next_slots = list(slots)
                    next_slots[axis] = sign
                    next_dp[(next_height, tuple(next_slots))] += count
            dp = next_dp
        return sum(dp.values())

    height_dp: dict[int, int] = {0: 1}
    for position, letters in enumerate(letters_by_position, start=1):
        next_height_dp: dict[int, int] = defaultdict(int)
        for height, count in height_dp.items():
            for letter in letters:
                next_height = height + DH[letter]
                if schedule_allows(next_height, position, schedule):
                    next_height_dp[next_height] += count
        height_dp = next_height_dp
    return sum(height_dp.values())


def l_formula(r: int) -> dict[str, int]:
    def suffix_count(length: int) -> int:
        return 1 if length == 0 else 2 * 3 ** (length - 1)

    corners = sum(3 ** k * suffix_count(r - 2 - k) for k in range(r - 1)) if r >= 2 else 0
    folds = sum(3 ** k * suffix_count(r - 3 - k) for k in range(r - 2)) if r >= 3 else 0
    return {"mono": 3 ** r, "G_each": corners, "D_each": folds,
            "L": 3 ** r + 4 * corners + 2 * folds}


def l_schedule_overlap(schedule: tuple[int, int, int, int, int],
                       coordinate_monotone: bool = False) -> tuple[int, dict[str, int]]:
    """Count each of the seven disjoint L languages with fixed defect slots."""
    if sum((schedule[0], schedule[2], schedule[4])) != R:
        raise ValueError("e120 intersections are defined at length 35")
    rows: dict[str, int] = {}
    rows["mono"] = language_schedule_count([M_LETTERS] * R, schedule, coordinate_monotone)

    for name, defect, predecessor, forbidden_successor in (
        ("G1", 3, 4, 2), ("G2", 5, 2, 4), ("G3", 3, 1, 2), ("G4", 5, 1, 4),
    ):
        count = 0
        for defect_position in range(1, R):
            language = [M_LETTERS] * R
            language[defect_position - 1] = (predecessor,)
            language[defect_position] = (defect,)
            if defect_position + 1 < R:
                language[defect_position + 1] = tuple(
                    letter for letter in M_LETTERS if letter != forbidden_successor
                )
            count += language_schedule_count(language, schedule, coordinate_monotone)
        rows[name] = count

    for name, middle in (("D1", 2), ("D2", 4)):
        count = 0
        for fold_position in range(2, R):
            language = [M_LETTERS] * R
            language[fold_position - 2] = (1,)
            language[fold_position - 1] = (middle,)
            language[fold_position] = (0,)
            if fold_position + 1 < R:
                language[fold_position + 1] = (2, 4)
            count += language_schedule_count(language, schedule, coordinate_monotone)
        rows[name] = count
    return sum(rows.values()), rows


# ---------------------------------------------------------------------------
# Prescribed two-defect languages.
# ---------------------------------------------------------------------------


def gadget_count(r: int) -> int:
    """Count either two-gadget language exactly by segment lengths."""
    if r < 7:
        return 0
    free_total = r - 6
    total = 0
    for first_length in range(free_total + 1):
        for middle_length in range(1, free_total - first_length + 1):
            last_length = free_total - first_length - middle_length
            middle_count = 2 * 3 ** (middle_length - 1)
            last_count = 1 if last_length == 0 else 2 * 3 ** (last_length - 1)
            total += 3 ** first_length * middle_count * last_count
    return total


def gadget_formula_35() -> int:
    n = R - 6
    return 2 * n * 3 ** (n - 1) + 2 * n * (n - 1) * 3 ** (n - 2)


def gadget_schedule_overlap(schedule: tuple[int, int, int, int, int], kind: str) -> int:
    """Exact finite-language H∩F or H∩Q count by free-segment lengths."""
    if kind == "fold":
        gadget = ((1,), (2,), (0,))       # -e1,+e2,+e1
        safe = (2, 4)                     # never immediately return left
    elif kind == "corner":
        gadget = ((2,), (4,), (3,))       # +e2,+e3,-e2
        safe = (1, 4)                     # never immediately return in +e2
    else:
        raise ValueError(kind)

    free_total = R - 6
    total = 0
    for first_length in range(free_total + 1):
        for middle_length in range(1, free_total - first_length + 1):
            last_length = free_total - first_length - middle_length
            language: list[tuple[int, ...]] = [M_LETTERS] * first_length
            language.extend(gadget)
            language.append(safe)
            language.extend([M_LETTERS] * (middle_length - 1))
            language.extend(gadget)
            if last_length:
                language.append(safe)
                language.extend([M_LETTERS] * (last_length - 1))
            if len(language) != R:
                raise AssertionError("gadget language length")
            total += language_schedule_count(language, schedule)
    return total


def one_defect_schedule_overlap(schedule: tuple[int, int, int, int, int], defect: int,
                                free: tuple[int, int]) -> int:
    """Exact H overlap for R_y/R_z, which are recorded as C-redundant."""
    total = 0
    for position in range(R):
        language = [free] * R
        language[position] = (defect,)
        total += language_schedule_count(language, schedule, coordinate_monotone=True)
    return total


# ---------------------------------------------------------------------------
# Full bounded split-search tables.
# ---------------------------------------------------------------------------


def two_block_table(blocks: list[dict[int, int]], tail: list[int]) -> list[dict[str, int]]:
    """All declared profile-feasible two-block schedules, not just the winner."""
    rows: list[dict[str, int]] = []
    for l1 in range(9, 14):
        for l2 in range(9, 14):
            tail_length = R - l1 - l2
            if not 9 <= tail_length <= 13:
                continue
            for q1, b1 in sorted(blocks[l1].items()):
                for q2, b2 in sorted(blocks[l2].items()):
                    rows.append({"l1": l1, "q1": q1, "l2": l2, "q2": q2, "t": tail_length,
                                 "B1": b1, "B2": b2, "T": tail[tail_length],
                                 "product": b1 * b2 * tail[tail_length]})
    return rows


def three_block_table(blocks: list[dict[int, int]], tail: list[int]) -> list[dict[str, int]]:
    """All declared profile-feasible asymmetric three-block schedules."""
    rows: list[dict[str, int]] = []
    for l1 in range(5, 10):
        for l2 in range(5, 10):
            for l3 in range(5, 10):
                tail_length = R - l1 - l2 - l3
                if not 8 <= tail_length <= 13:
                    continue
                for q1, b1 in sorted(blocks[l1].items()):
                    for q2, b2 in sorted(blocks[l2].items()):
                        for q3, b3 in sorted(blocks[l3].items()):
                            rows.append({
                                "l1": l1, "q1": q1, "l2": l2, "q2": q2,
                                "l3": l3, "q3": q3, "t": tail_length,
                                "B1": b1, "B2": b2, "B3": b3, "T": tail[tail_length],
                                "product": b1 * b2 * b3 * tail[tail_length],
                            })
    return rows


# ---------------------------------------------------------------------------
# Small finite check of the first-backtrack chain implementation.
# ---------------------------------------------------------------------------


def saw_series(max_depth: int, forbidden: tuple[int, ...] = ()) -> list[int]:
    counts = [0] * (max_depth + 1)
    origin = pack(0, 0, 0)
    forbidden_set = frozenset(forbidden)
    visited = {origin}

    def visit(point: int, depth: int) -> None:
        counts[depth] += 1
        if depth == max_depth:
            return
        for delta in PACKED_STEPS:
            nxt = point + delta
            if nxt in forbidden_set or nxt in visited:
                continue
            visited.add(nxt)
            visit(nxt, depth + 1)
            visited.remove(nxt)

    visit(origin, 0)
    return counts


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
        context.prec = 100
        return format(Decimal(value).quantize(Decimal(1).scaleb(-REPORT_PLACES), rounding=ROUND_FLOOR), "f")


def main() -> None:
    started = time.monotonic()
    timings: dict[str, float] = {}

    source_hash = hashlib.sha256(SAW_SOURCE.read_bytes()).hexdigest()
    check("external_c36_source_hash", source_hash == C36_SHA256,
          "[EXTERNAL] cached Schram--Barkema--Bisseling Table-I PDF hash")

    t0 = time.monotonic()
    tail, blocks = positive_height_profile(13)
    timings["positive_height_profile"] = time.monotonic() - t0
    check("T13", tail[13] == 142016661,
          "[COMPUTATION] exact visited-set profile gives T(13)=142016661")
    check("B11_profile", blocks[11] == {3: 5832, 5: 411156, 7: 729000, 9: 551124, 11: 177147},
          "[COMPUTATION] all five nonzero B(11,q) entries from one record-height DFS")

    t0 = time.monotonic()
    two_rows = two_block_table(blocks, tail)
    three_rows = three_block_table(blocks, tail)
    timings["split_table_assembly"] = time.monotonic() - t0
    two_max = max(row["product"] for row in two_rows)
    three_max = max(row["product"] for row in three_rows)
    two_winners = [row for row in two_rows if row["product"] == two_max]
    three_winners = [row for row in three_rows if row["product"] == three_max]
    check("two_block_full_search", len(two_rows) == 392 and len(two_winners) == 1
          and two_winners[0] == {"l1": 11, "q1": 7, "l2": 11, "q2": 7, "t": 13,
                                  "B1": 729000, "B2": 729000, "T": 142016661,
                                  "product": 75473476338501000000},
          "[COMPUTATION] all 392 declared two-block schedules; unique product maximum is wave-11 H")
    check("three_block_full_search", len(three_rows) == 1651 and len(three_winners) == 6
          and three_max == 17606452560245513280,
          "[COMPUTATION] all 1651 declared asymmetric three-block schedules; maximum is below H")

    l_data = l_formula(R)
    l35 = l_data["L"]
    c35 = cm_automaton_count(R)
    c_cap_l = 3 ** R + 4 * (R - 1) * 2 ** (R - 2)
    check("L35", l35 == 644233352324156721,
          "[THEOREM] re-derived seven-family L formula at r=35")
    check("C35", c35 == 333543290395947705,
          "[COMPUTATION] 27-state signed-coordinate automaton counts C_35")
    check("C_cap_L", c_cap_l == 50032713330104219,
          "[THEOREM] C∩L = 3^35 + 4·34·2^33")

    # The 25 schedules H_{q1,q2} are pairwise disjoint: their values at the
    # fixed times 11 and 22 are (q1,q1+q2), respectively.
    t0 = time.monotonic()
    grid_rows: list[dict[str, object]] = []
    mobius_tables: dict[str, list[dict[str, object]]] = {}
    for q1 in sorted(blocks[11]):
        for q2 in sorted(blocks[11]):
            schedule = (11, q1, 11, q2, 13)
            h_count = blocks[11][q1] * blocks[11][q2] * tail[13]
            c_overlap, table = cm_overlap_mobius(schedule)
            l_overlap, l_parts = l_schedule_overlap(schedule, coordinate_monotone=False)
            cl_overlap, cl_parts = l_schedule_overlap(schedule, coordinate_monotone=True)
            key = f"{q1},{q2}"
            mobius_tables[key] = table
            grid_rows.append({
                "q1": q1,
                "q2": q2,
                "schedule": list(schedule),
                "H": h_count,
                "C_cap_H": c_overlap,
                "L_cap_H": l_overlap,
                "C_cap_L_cap_H": cl_overlap,
                "L_parts": l_parts,
                "CL_parts": cl_parts,
            })
    timings["grid_intersections"] = time.monotonic() - t0
    grid_c_sum = sum(int(row["C_cap_H"]) for row in grid_rows)
    grid_l_sum = sum(int(row["L_cap_H"]) for row in grid_rows)
    grid_cl_sum = sum(int(row["C_cap_L_cap_H"]) for row in grid_rows)
    check("height_grid_pairwise_disjoint", len(grid_rows) == 25,
          "[LEMMA] distinct (q1,q2) prescribe distinct heights at steps 11 or 22")
    check("height_grid_overlap_totals", grid_c_sum == 51027099125602155
          and grid_l_sum == 491050350045737865
          and grid_cl_sum == 50032404092458907,
          "[COMPUTATION] exact Möbius/regular-language pairwise and triple intersection totals")

    # Candidate R_y/R_z are recorded rather than added because each is already
    # a coordinate-monotone C subfamily.
    r_one_defect = R * 2 ** (R - 1)
    ry_h_rows: dict[str, int] = {}
    rz_h_rows: dict[str, int] = {}
    for row in grid_rows:
        q1, q2 = int(row["q1"]), int(row["q2"])
        schedule = (11, q1, 11, q2, 13)
        key = f"{q1},{q2}"
        ry_h_rows[key] = one_defect_schedule_overlap(schedule, 3, (1, 4))
        rz_h_rows[key] = one_defect_schedule_overlap(schedule, 5, (1, 2))
    check("one_defect_candidate_audit", r_one_defect == 601295421440
          and sum(ry_h_rows.values()) == sum(rz_h_rows.values()) == 429496729600,
          "[THEOREM] R_y,R_z have size 35·2^34, lie in C, and cannot enlarge the union")

    # Two new disjoint languages.  Each contains two non-M defects, hence is
    # disjoint from L; each reverses a coordinate sign, hence is disjoint from
    # C; their letter-0 content distinguishes fold from corner.
    t0 = time.monotonic()
    gadget35 = gadget_count(R)
    fold_h_rows: dict[str, int] = {}
    corner_h_rows: dict[str, int] = {}
    for row in grid_rows:
        q1, q2 = int(row["q1"]), int(row["q2"])
        schedule = (11, q1, 11, q2, 13)
        key = f"{q1},{q2}"
        fold_h_rows[key] = gadget_schedule_overlap(schedule, "fold")
        corner_h_rows[key] = gadget_schedule_overlap(schedule, "corner")
    timings["two_defect_intersections"] = time.monotonic() - t0
    fold_h_sum = sum(fold_h_rows.values())
    corner_h_sum = sum(corner_h_rows.values())
    check("two_defect_language_count", gadget35 == gadget_formula_35() == 13710824278006626,
          "[THEOREM] each prescribed two-defect language has 2N·3^(N-1)+2N(N-1)·3^(N-2), N=29")
    check("two_defect_intersections", fold_h_sum == corner_h_sum == 7732355849776818,
          "[COMPUTATION] exact finite-language intersections with all 25 H schedules")

    # One safe y-fold is the pre-existing D1 language.  Its exact audit shows
    # why the requested one-fold candidate contributes no new set mass.
    d1_h_rows = {f"{row['q1']},{row['q2']}": int(dict(row["L_parts"])["D1"]) for row in grid_rows}
    check("one_fold_candidate_redundant", sum(d1_h_rows.values()) == 31501343210481297,
          "[THEOREM] the one-fold candidate is D1⊂L; all its H intersections are retained")

    grid_net = sum(int(row["H"]) - int(row["C_cap_H"]) - int(row["L_cap_H"])
                   + int(row["C_cap_L_cap_H"]) for row in grid_rows)
    fold_net = gadget35 - fold_h_sum
    corner_net = gadget35 - corner_h_sum
    union = l35 + c35 - c_cap_l + grid_net + fold_net + corner_net
    m_value = C36 - union
    check("union_a35", union == 499330428831189067251,
          "[THEOREM] exact inclusion--exclusion gives the enlarged finite a_35 family")
    check("M_value", m_value == 2940871525905870537493419,
          "[THEOREM] M=c_36-a_35_lower is strictly below the wave-11 M")
    check("strict_improvement", union > 76401062626946721319
          and m_value < 2941294455272074779839351,
          "[THEOREM] the enlarged union strictly improves the prior endpoint")

    # A finite brute-force implementation check of the all-size first-backtrack
    # injection.  The all-size proof is written in proofs/saw_union2.md.
    t0 = time.monotonic()
    c_small = saw_series(8)
    a_small = saw_series(8, (pack(1, 0, 0),))
    chain_rows: list[dict[str, int | bool]] = []
    for total in range(2, 9):
        for n in range(1, total):
            m = total - n
            rhs = c_small[m] * (c_small[n] - a_small[n - 1])
            chain_rows.append({"m": m, "n": n, "lhs": c_small[total], "rhs": rhs,
                               "passed": c_small[total] <= rhs})
    timings["small_chain_check"] = time.monotonic() - t0
    check("first_backtrack_small", all(bool(row["passed"]) for row in chain_rows),
          "[COMPUTATION] c_(m+n)<=c_m(c_n-a_(n-1)) for all m+n<=8")

    lower, upper = atanh_root_interval(m_value)
    floor = floor40(lower)
    check("endpoint_floor", floor == "0.2122129322754723621039731646196408804906",
          "[COMPUTATION] directed 120-dps floor of atanh(M^(-1/36))")
    with localcontext() as context:
        context.prec = 180
        check("endpoint_rounding_semantics", Decimal(floor) <= Decimal(lower)
              and Decimal(upper) < Decimal(floor) + Decimal(1).scaleb(-REPORT_PLACES),
              "[COMPUTATION] floor<=interval lower<=true<=upper<floor+10^-40")

    elapsed = time.monotonic() - started
    payload = {
        "provenance": {
            "script": SCRIPT,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": ".venv/bin/python",
            "arithmetic": "[COMPUTATION] Python integers for all finite counts; mpmath.iv directed interval arithmetic at 120 dps only for the endpoint.",
            "mpmath_iv_dps": DPS,
            "elapsed_seconds": elapsed,
            "timings": timings,
        },
        "data": {
            "classification": "[THEOREM] a_35 >= 499330428831189067251 by an exact finite union of L, C, 25 pairwise-disjoint height schedules, and two disjoint prescribed two-defect languages; [UNRESOLVED] this is a finite n=35 certificate, not an all-size lower-family theorem.",
            "notation": {
                "c_r": "rooted r-step cubic-lattice SAW count",
                "a_r": "rooted r-step SAW count avoiding e_1=(1,0,0)",
                "h": "h(x,y,z)=-x+y+z; letter increments (-1,+1,+1,-1,+1,-1)",
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
                "method": "visited-set DFS of positive-height SAWs; a strict-record terminal q is exactly B(n,q)",
                "tail_T_0_to_13": tail,
                "blocks_B_n_q": {str(n): {str(q): count for q, count in sorted(blocks[n].items())}
                                  for n in range(1, 14)},
            },
            "split_search": {
                "status": "[COMPUTATION]",
                "two_block_domain": "all l1,l2,t in [9,13] with l1+l2+t=35 and every nonzero B(li,qi) from the stored profile",
                "three_block_domain": "all l1,l2,l3 in [5,9], t in [8,13], l1+l2+l3+t=35 and every nonzero B(li,qi) from the stored profile",
                "two_block_table": two_rows,
                "three_block_table": three_rows,
                "two_block_winner": two_winners[0],
                "three_block_max_product": three_max,
                "three_block_winners": three_winners,
                "selection": "[COMPUTATION] the old (11,7,11,7,13) product is the unique two-block maximum in the declared exact profile search; all 25 (11,q1,11,q2,13) schedules are nevertheless used because they are pairwise height-disjoint.",
            },
            "base_families": {
                "L": {"status": "[THEOREM]", "count": l35, "formula": "3^r+4G_r+2D_r"},
                "C": {"status": "[COMPUTATION]", "count": c35,
                      "automaton": "27 sign-slot states; first +e1 forbidden; independently recomputed by test per used signed pattern"},
                "C_cap_L": {"status": "[THEOREM]", "count": c_cap_l,
                            "formula": "3^35+4·34·2^33"},
            },
            "height_grid": {
                "status": "[THEOREM]",
                "definition": "H_(q1,q2): B(11,q1), translated B(11,q2), then T(13), for q1,q2 in {3,5,7,9,11}",
                "injection": "[THEOREM] strict disjoint height bands make every concatenation self-avoiding and recoverable at steps 11,22; positive height avoids e_1.",
                "disjointness": "[LEMMA] different q1 or q2 prescribe different heights at step 11 or 22, so the 25 schedules are pairwise disjoint.",
                "rows": grid_rows,
                "cm_mobius_tables": mobius_tables,
                "overlap_method": "[THEOREM] signed-orthant Möbius inversion for C∩H; seven regular-language height DPs for L∩H and C∩L∩H.",
                "totals": {"H_sum": sum(int(row["H"]) for row in grid_rows), "C_cap_H_sum": grid_c_sum,
                           "L_cap_H_sum": grid_l_sum, "C_cap_L_cap_H_sum": grid_cl_sum},
            },
            "candidate_audit": {
                "one_defect_coordinate_monotone": {
                    "status": "[THEOREM] redundant",
                    "R_y_count": r_one_defect,
                    "R_z_count": r_one_defect,
                    "proof": "each uses exactly one -e2 or -e3 and one sign per coordinate, hence is a C subset; R_y∩L=R_z∩L=34·2^33 and R_y∩R_z=empty.",
                    "R_y_cap_H_grid": ry_h_rows,
                    "R_z_cap_H_grid": rz_h_rows,
                },
                "one_fold": {
                    "status": "[THEOREM] redundant",
                    "definition": "u(-e1,+e2,+e1)w with u,w in M* and w empty or beginning +e2/+e3",
                    "identification": "this is exactly the existing D1 language, so it is contained in L and adds no union mass",
                    "D1_cap_H_grid": d1_h_rows,
                },
            },
            "two_defect_families": {
                "fold": {
                    "status": "[THEOREM]",
                    "definition": "u F_y v F_y w; F_y=(-e1,+e2,+e1); u∈M*, v∈M+ starts +e2/+e3, w empty or starts +e2/+e3",
                    "count": gadget35,
                    "formula": "N=29: 2N·3^(N-1)+2N(N-1)·3^(N-2)",
                    "self_avoidance": "[LEMMA] each gadget is a one-square detour; the required first post-gadget advance prevents its unique possible return, and monotonicity thereafter prevents every other collision.",
                    "e1_avoidance": "[LEMMA] the x coordinate is never positive.",
                    "disjointness": "[THEOREM] two +e1 defect letters exclude L; both x signs exclude C; corner family contains no +e1 and is disjoint.",
                    "cap_H_grid": fold_h_rows,
                    "cap_H_grid_sum": fold_h_sum,
                },
                "corner": {
                    "status": "[THEOREM]",
                    "definition": "u G_y v G_y w; G_y=(+e2,+e3,-e2); u∈M*, v∈M+ starts -e1/+e3, w empty or starts -e1/+e3",
                    "count": gadget35,
                    "formula": "N=29: 2N·3^(N-1)+2N(N-1)·3^(N-2)",
                    "self_avoidance": "[LEMMA] each gadget is a one-square detour; the required first post-gadget advance prevents its unique possible return, and monotonicity thereafter prevents every other collision.",
                    "e1_avoidance": "[LEMMA] no +e1 letter occurs.",
                    "disjointness": "[THEOREM] two -e2 defect letters exclude L; both y signs exclude C; fold family contains +e1 and is disjoint.",
                    "cap_H_grid": corner_h_rows,
                    "cap_H_grid_sum": corner_h_sum,
                },
            },
            "union_certificate": {
                "status": "[THEOREM]",
                "expression": "|L|+|C|-|L∩C|+sum_s(|H_s|-|H_s∩C|-|H_s∩L|+|H_s∩C∩L|)+|F|-sum_s|F∩H_s|+|Q|-sum_s|Q∩H_s|",
                "a35_lower": union,
                "components": {"L": l35, "C": c35, "C_cap_L": c_cap_l,
                               "height_grid_net": grid_net, "fold_net": fold_net,
                               "corner_net": corner_net},
                "scope": "[UNRESOLVED] a finite n=35 subfamily statement only; no all-r strengthening is claimed.",
            },
            "chain": {
                "status": "[THEOREM]",
                "statement": "for all m,n>=1, c_(m+n)<=c_m(c_n-a_(n-1)); hence mu^n<=c_n-a_(n-1)",
                "n36_instance": "mu^36<=c_36-a_35_lower",
                "M": m_value,
                "small_exact_rows": chain_rows,
            },
            "endpoint": {
                "status": "[THEOREM] via the SAW-domination proof in proofs/kc_bounds.md §A",
                "statement": "K_c>=atanh(M^(-1/36))",
                "interval_120dps": [lower, upper],
                "floor_40": floor,
            },
            "resource_walls": [],
            "unresolved": [
                "[UNRESOLVED] the bounded 392/1651-row split searches are exact over their stated profile-feasible domains, not global optimizations over arbitrary 35-step schedules.",
                "[UNRESOLVED] no all-size lower-family recursion beyond the existing L_r theorem is established.",
                "[UNRESOLVED] broader one-fold candidates are not claimed: the explicitly audited y-fold is D1 subset L.",
            ],
        },
        "checks": checks,
    }
    RESULT.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    print(f"\nwrote {RESULT}")
    print(f"elapsed {elapsed:.1f}s")
    print(f"a_35 lower bound: {union}")
    print(f"M: {m_value}")
    print(f"K_c lower floor (40): {floor}")


if __name__ == "__main__":
    try:
        main()
    except AssertionError as exc:
        print(f"CERTIFICATE FAILURE: {exc}")
        raise
