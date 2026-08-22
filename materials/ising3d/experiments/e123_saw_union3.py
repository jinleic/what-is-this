#!/usr/bin/env python3
"""Exact variable-cut and three-block SAW-union certificate at n=36.

[THEOREM] This producer extends the already certified finite a_35 union by
four variable-cut two-block height slabs, twelve three-block height slabs, and
two new three-gadget defect languages.  All finite counts are Python integers.
The only noninteger operation is the recorded directed interval enclosure of
the resulting Ising lower endpoint.

Run:
    PYTHONPATH=src .venv/bin/python experiments/e123_saw_union3.py
    PYTHONPATH=src .venv/bin/python tests/test_saw_union3.py
"""
from __future__ import annotations

from collections import defaultdict
from datetime import datetime, timezone
from decimal import Decimal, ROUND_FLOOR, localcontext
from itertools import product
import hashlib
import json
from math import prod
from pathlib import Path
import resource
import sys
import time

import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = "experiments/e123_saw_union3.py"
RESULT = ROOT / "results" / "bounds" / "saw_union3.json"
SAW_SOURCE = ROOT / "sources" / "fulltext" / "schram_barkema_bisseling2011.pdf"

C36 = 2941370856334701726560670
C36_SHA256 = "898d8333e8d70acd279e29f226ff108550b132463aa662f86df590b0ad72cb12"
DPS = 120
REPORT_PLACES = 40
R = 35
M = (1, 2, 4)  # {-e1,+e2,+e3}; all have +1 height increment.
STEPS = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))
DH = (-1, 1, 1, -1, 1, -1)
OFF = 40

OLD_GRID_Q = (3, 5, 7, 9, 11)
VARIABLE_SCHEDULES = (
    ("V5", (13, 5, 11, 3, 11)),
    ("V7", (13, 7, 11, 3, 11)),
    ("V9", (13, 9, 11, 3, 11)),
    ("V11", (13, 11, 11, 3, 11)),
)
THREE_SCHEDULES = tuple(
    (f"J{q1}_{q3}", (7, q1, 7, 3, 9, q3, 12))
    for q1 in (3, 5, 7) for q3 in (3, 5, 7, 9)
)

checks: list[dict[str, object]] = []


def check(name: str, passed: bool, detail: str) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})
    if not passed:
        print(f"FAIL {name}: {detail}")
        raise AssertionError(name)
    print(f"PASS {name}: {detail}")


# ---------------------------------------------------------------------------
# Packed-coordinate raw positive-height SAW enumeration.
# ---------------------------------------------------------------------------


def pack(x: int, y: int, z: int) -> int:
    return ((x + OFF) << 14) | ((y + OFF) << 7) | (z + OFF)


def pack_delta(dx: int, dy: int, dz: int) -> int:
    return (dx << 14) + (dy << 7) + dz


PACKED_STEPS = tuple(pack_delta(*step) for step in STEPS)


def positive_height_profile(max_depth: int) -> tuple[list[int], list[dict[int, int]]]:
    """Return T(n) and strict-record B(n,q) through ``max_depth``.

    [LEMMA] A positive-height terminal at q is a B(n,q) block precisely if q
    is a strict record height.  The recursive packed DFS visits each raw SAW
    once and is intentionally unlike the test's tuple iterative traversal.
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
# Height schedules, direct separators, and exact intersection DPs.
# ---------------------------------------------------------------------------


def allowed_heights(schedule: tuple[int, ...], position: int) -> set[int]:
    """Return the finite set of physically possible allowed heights at a step."""
    if len(schedule) == 5:
        l1, q1, l2, q2, _tail = schedule
        if position < l1:
            return set(range(1, q1))
        if position == l1:
            return {q1}
        if position < l1 + l2:
            return set(range(q1 + 1, q1 + q2))
        if position == l1 + l2:
            return {q1 + q2}
        return set(range(q1 + q2 + 1, R + 1))
    l1, q1, l2, q2, l3, q3, _tail = schedule
    cut2 = l1 + l2
    cut3 = cut2 + l3
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


def schedule_allows(schedule: tuple[int, ...], height: int, position: int) -> bool:
    return height in allowed_heights(schedule, position)


def first_separator(left: tuple[int, ...], right: tuple[int, ...]) -> dict[str, object] | None:
    """An empty allowed-height intersection is a direct H-left/H-right proof."""
    for position in range(1, R + 1):
        left_heights = allowed_heights(left, position)
        right_heights = allowed_heights(right, position)
        if not left_heights.intersection(right_heights):
            return {
                "step": position,
                "left_allowed": sorted(left_heights),
                "right_allowed": sorted(right_heights),
            }
    return None


def schedule_components(schedule: tuple[int, ...]) -> tuple[tuple[tuple[int, int], ...], int]:
    if len(schedule) == 5:
        l1, q1, l2, q2, tail_length = schedule
        return ((l1, q1), (l2, q2)), tail_length
    l1, q1, l2, q2, l3, q3, tail_length = schedule
    return ((l1, q1), (l2, q2), (l3, q3)), tail_length


def block_word_dp(length: int, rise: int, letters: tuple[int, ...]) -> int:
    """Exact height DP for a coordinate-monotone strict-record block."""
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


def tail_word_dp(length: int, letters: tuple[int, ...]) -> int:
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


def orthant_letters(code: tuple[int | None, int | None, int | None]) -> tuple[int, ...]:
    return tuple(axis * 2 + (0 if sign == 0 else 1)
                 for axis, sign in enumerate(code) if sign is not None)


def cm_overlap_mobius(schedule: tuple[int, ...]) -> tuple[int, list[dict[str, object]]]:
    """[THEOREM] Count C∩H by exact signed-orthant Möbius inversion."""
    block_components, tail_length = schedule_components(schedule)
    rows: list[dict[str, object]] = []
    for code in product((None, 0, 1), repeat=3):
        letters = orthant_letters(code)
        if not letters:
            continue
        empty_coordinates = sum(sign is None for sign in code)
        block_counts = [block_word_dp(length, rise, letters) for length, rise in block_components]
        tail_count = tail_word_dp(tail_length, letters)
        term = (-1) ** empty_coordinates * prod(block_counts) * tail_count
        rows.append({
            "code": list(code),
            "letters": list(letters),
            "empty_coordinates": empty_coordinates,
            "sign": (-1) ** empty_coordinates,
            "block_counts": block_counts,
            "tail": tail_count,
            "term": term,
        })
    return sum(int(row["term"]) for row in rows), rows


def words_in_schedule(letters_by_position: list[tuple[int, ...]], schedule: tuple[int, ...],
                      coordinate_monotone: bool = False) -> int:
    """Exact finite language DP; sign slots are retained only for C∩L∩H."""
    if coordinate_monotone:
        dp: dict[tuple[int, tuple[int, int, int]], int] = {(0, (0, 0, 0)): 1}
        for position, letters in enumerate(letters_by_position, start=1):
            next_dp: dict[tuple[int, tuple[int, int, int]], int] = defaultdict(int)
            for (height, slots), count in dp.items():
                for letter in letters:
                    next_height = height + DH[letter]
                    if not schedule_allows(schedule, next_height, position):
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

    dp_height: dict[int, int] = {0: 1}
    for position, letters in enumerate(letters_by_position, start=1):
        next_dp_height: dict[int, int] = defaultdict(int)
        for height, count in dp_height.items():
            for letter in letters:
                next_height = height + DH[letter]
                if schedule_allows(schedule, next_height, position):
                    next_dp_height[next_height] += count
        dp_height = next_dp_height
    return sum(dp_height.values())


# ---------------------------------------------------------------------------
# Existing L/C family intersections.
# ---------------------------------------------------------------------------


def l_formula(r: int) -> dict[str, int]:
    def suffix_count(length: int) -> int:
        return 1 if length == 0 else 2 * 3 ** (length - 1)

    corners = sum(3 ** k * suffix_count(r - 2 - k) for k in range(r - 1)) if r >= 2 else 0
    folds = sum(3 ** k * suffix_count(r - 3 - k) for k in range(r - 2)) if r >= 3 else 0
    return {"mono": 3 ** r, "G_each": corners, "D_each": folds,
            "L": 3 ** r + 4 * corners + 2 * folds}


def cm_automaton_count(length: int) -> int:
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


def l_schedule_overlap(schedule: tuple[int, ...], coordinate_monotone: bool = False) -> tuple[int, dict[str, int]]:
    """Count all seven disjoint regular L languages against one height schedule."""
    rows: dict[str, int] = {}
    rows["mono"] = words_in_schedule([M] * R, schedule, coordinate_monotone)
    for name, defect, predecessor, forbidden_successor in (
        ("G1", 3, 4, 2), ("G2", 5, 2, 4), ("G3", 3, 1, 2), ("G4", 5, 1, 4),
    ):
        count = 0
        for defect_position in range(1, R):
            language = [M] * R
            language[defect_position - 1] = (predecessor,)
            language[defect_position] = (defect,)
            if defect_position + 1 < R:
                language[defect_position + 1] = tuple(letter for letter in M if letter != forbidden_successor)
            count += words_in_schedule(language, schedule, coordinate_monotone)
        rows[name] = count
    for name, middle in (("D1", 2), ("D2", 4)):
        count = 0
        for fold_position in range(2, R):
            language = [M] * R
            language[fold_position - 2] = (1,)
            language[fold_position - 1] = (middle,)
            language[fold_position] = (0,)
            if fold_position + 1 < R:
                language[fold_position + 1] = (2, 4)
            count += words_in_schedule(language, schedule, coordinate_monotone)
        rows[name] = count
    return sum(rows.values()), rows


# ---------------------------------------------------------------------------
# Prescribed k-gadget languages and their schedule intersections.
# ---------------------------------------------------------------------------


def gadget_spec(kind: str) -> tuple[tuple[int, int, int], tuple[int, int]]:
    if kind == "fold":
        return (1, 2, 0), (2, 4)
    if kind == "corner":
        return (2, 4, 3), (1, 4)
    raise ValueError(kind)


def free_segment_lengths(total: int, gadgets: int):
    """Yield n0,n1,...,n_(g-1),ng with internal segments nonempty."""
    for first in range(total + 1):
        def append_internal(prefix: tuple[int, ...], remaining: int, pending: int):
            if pending == 0:
                yield prefix + (remaining,)
                return
            for length in range(1, remaining - (pending - 1) + 1):
                yield from append_internal(prefix + (length,), remaining - length, pending - 1)
        yield from append_internal((first,), total - first, gadgets - 1)


def gadget_count(r: int, gadgets: int) -> int:
    free_total = r - 3 * gadgets
    if free_total < gadgets - 1:
        return 0
    total = 0
    for lengths in free_segment_lengths(free_total, gadgets):
        first_length, *rest_lengths = lengths
        factor = 3 ** first_length
        for length in rest_lengths[:-1]:
            factor *= 2 * 3 ** (length - 1)
        final_length = rest_lengths[-1]
        factor *= 1 if final_length == 0 else 2 * 3 ** (final_length - 1)
        total += factor
    return total


def gadget_schedule_overlap(schedule: tuple[int, ...], kind: str, gadgets: int) -> int:
    """Exact finite word DP for a k-gadget language ∩ one height schedule."""
    free_total = R - 3 * gadgets
    if free_total < gadgets - 1:
        return 0
    gadget, safe = gadget_spec(kind)
    singleton_gadget = [(letter,) for letter in gadget]
    total = 0
    for lengths in free_segment_lengths(free_total, gadgets):
        first_length, *rest_lengths = lengths
        internal_lengths = rest_lengths[:-1]
        final_length = rest_lengths[-1]
        language: list[tuple[int, ...]] = [M] * first_length
        for index in range(gadgets):
            language.extend(singleton_gadget)
            if index < gadgets - 1:
                length = internal_lengths[index]
                language.append(safe)
                language.extend([M] * (length - 1))
            elif final_length:
                language.append(safe)
                language.extend([M] * (final_length - 1))
        if len(language) != R:
            raise AssertionError("gadget language length")
        total += words_in_schedule(language, schedule)
    return total


# ---------------------------------------------------------------------------
# Exact quantities and artifact assembly.
# ---------------------------------------------------------------------------


def height_schedule_count(schedule: tuple[int, ...], blocks: list[dict[int, int]], tail: list[int]) -> int:
    components, tail_length = schedule_components(schedule)
    return prod(blocks[length][rise] for length, rise in components) * tail[tail_length]


def build_schedule_row(identifier: str, schedule: tuple[int, ...], category: str,
                       blocks: list[dict[int, int]], tail: list[int],
                       mobius_tables: dict[str, list[dict[str, object]]]) -> dict[str, object]:
    c_overlap, mobius_table = cm_overlap_mobius(schedule)
    l_overlap, _ = l_schedule_overlap(schedule)
    cl_overlap, _ = l_schedule_overlap(schedule, coordinate_monotone=True)
    mobius_tables[identifier] = mobius_table
    return {
        "id": identifier,
        "schedule": list(schedule),
        "H": height_schedule_count(schedule, blocks, tail),
        "C_cap_H": c_overlap,
        "L_cap_H": l_overlap,
        "C_cap_L_cap_H": cl_overlap,
        "F2_cap_H": gadget_schedule_overlap(schedule, "fold", 2),
        "Q2_cap_H": gadget_schedule_overlap(schedule, "corner", 2),
        **({"category": category} if category else {}),
    }


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
        context.prec = 180
        return format(Decimal(value).quantize(Decimal(1).scaleb(-REPORT_PLACES), rounding=ROUND_FLOOR), "f")


def peak_rss_bytes() -> int:
    raw = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(raw if sys.platform == "darwin" else raw * 1024)


def main() -> None:
    started = time.monotonic()
    timings: dict[str, float] = {}

    source_hash = hashlib.sha256(SAW_SOURCE.read_bytes()).hexdigest()
    check("external_c36_source_hash", source_hash == C36_SHA256,
          "[EXTERNAL] cached Schram--Barkema--Bisseling Table-I PDF hash")

    t0 = time.monotonic()
    tail, blocks = positive_height_profile(13)
    timings["positive_height_profile"] = time.monotonic() - t0
    check("positive_height_profile", tail == [1, 3, 9, 45, 171, 837, 3411, 16425, 69525, 331947,
                                                 1436643, 6827355, 29971017, 142016661]
          and blocks[11] == {3: 5832, 5: 411156, 7: 729000, 9: 551124, 11: 177147},
          "[COMPUTATION] packed raw DFS reproduces all selected T/B values")

    l35 = l_formula(R)["L"]
    c35 = cm_automaton_count(R)
    c_cap_l = 3 ** R + 4 * (R - 1) * 2 ** (R - 2)
    check("base_L_C", l35 == 644233352324156721 and c35 == 333543290395947705
          and c_cap_l == 50032713330104219,
          "[THEOREM] existing L/C counts and C∩L formula at r=35")

    # Recompute the complete old base, retaining every term needed when the new
    # families are unioned relative to it.
    t0 = time.monotonic()
    mobius_tables: dict[str, list[dict[str, object]]] = {}
    old_rows: list[dict[str, object]] = []
    for q1 in OLD_GRID_Q:
        for q2 in OLD_GRID_Q:
            old_rows.append(build_schedule_row(f"H{q1}_{q2}", (11, q1, 11, q2, 13), "",
                                               blocks, tail, mobius_tables))
    timings["old_grid_terms"] = time.monotonic() - t0
    old_grid_c = sum(int(row["C_cap_H"]) for row in old_rows)
    old_grid_l = sum(int(row["L_cap_H"]) for row in old_rows)
    old_grid_cl = sum(int(row["C_cap_L_cap_H"]) for row in old_rows)
    old_f2 = sum(int(row["F2_cap_H"]) for row in old_rows)
    old_q2 = sum(int(row["Q2_cap_H"]) for row in old_rows)
    old_grid_net = sum(int(row["H"]) - int(row["C_cap_H"]) - int(row["L_cap_H"])
                       + int(row["C_cap_L_cap_H"]) for row in old_rows)
    two_gadget = gadget_count(R, 2)
    old_base = l35 + c35 - c_cap_l + old_grid_net + 2 * (two_gadget - old_f2)
    check("old_union_terms", len(old_rows) == 25 and old_grid_c == 51027099125602155
          and old_grid_l == 491050350045737865 and old_grid_cl == 50032404092458907
          and two_gadget == 13710824278006626 and old_f2 == old_q2 == 7732355849776818
          and old_base == 499330428831189067251,
          "[COMPUTATION] exact inherited H/C/L/F2/Q2 inclusion--exclusion base")

    # All 41 height schedules are pairwise separated at a literal time slice.
    # This provides each zero H_i∩H_j term; all higher grid intersections are
    # consequently zero as subsets of one of these listed pairs.
    schedule_specs = [(str(row["id"]), tuple(row["schedule"])) for row in old_rows]
    schedule_specs += list(VARIABLE_SCHEDULES) + list(THREE_SCHEDULES)
    separator_rows: list[dict[str, object]] = []
    for left_index, (left_id, left_schedule) in enumerate(schedule_specs):
        for right_id, right_schedule in schedule_specs[left_index + 1:]:
            separator = first_separator(left_schedule, right_schedule)
            if separator is None:
                raise AssertionError(f"unseparated schedules {left_id},{right_id}")
            separator_rows.append({"left": left_id, "right": right_id, **separator})
    check("all_schedule_pairs_disjoint", len(schedule_specs) == 41 and len(separator_rows) == 820,
          "[THEOREM] every pair has a recorded disjoint allowed-height slice")

    t0 = time.monotonic()
    variable_rows = [build_schedule_row(identifier, schedule, "variable_cut", blocks, tail, mobius_tables)
                     for identifier, schedule in VARIABLE_SCHEDULES]
    three_rows = [build_schedule_row(identifier, schedule, "three_block", blocks, tail, mobius_tables)
                  for identifier, schedule in THREE_SCHEDULES]
    timings["new_grid_terms"] = time.monotonic() - t0
    new_rows = variable_rows + three_rows

    def row_net(row: dict[str, object]) -> int:
        return (int(row["H"]) - int(row["C_cap_H"]) - int(row["L_cap_H"])
                + int(row["C_cap_L_cap_H"]) - int(row["F2_cap_H"]) - int(row["Q2_cap_H"]))

    variable_h = sum(int(row["H"]) for row in variable_rows)
    variable_c = sum(int(row["C_cap_H"]) for row in variable_rows)
    three_h = sum(int(row["H"]) for row in three_rows)
    three_c = sum(int(row["C_cap_H"]) for row in three_rows)
    variable_net = sum(row_net(row) for row in variable_rows)
    three_net = sum(row_net(row) for row in three_rows)
    check("variable_cut_grid", variable_h == 1366157335427744400 and variable_c == 4861170126720
          and all(int(row[key]) == 0 for row in variable_rows
                  for key in ("L_cap_H", "C_cap_L_cap_H", "F2_cap_H", "Q2_cap_H")),
          "[COMPUTATION] all four variable-cut schedule pair/triple terms")
    check("three_block_grid", three_h == 7141848962632019280 and three_c == 18012995095392
          and all(int(row[key]) == 0 for row in three_rows
                  for key in ("L_cap_H", "C_cap_L_cap_H", "F2_cap_H", "Q2_cap_H")),
          "[COMPUTATION] all twelve three-block schedule pair/triple terms")
    check("new_grid_net", variable_net + three_net == 8507983423894541568,
          "[THEOREM] exact new-grid contribution relative to the full old union")

    # F3 and Q3 are separately disjoint from L,C,F2,Q2 by their defect/sign
    # content.  Retain their intersections with every one of the 41 mutually
    # disjoint height schedules, including the sixteen exact zeros in new rows.
    t0 = time.monotonic()
    all_grid_rows = old_rows + new_rows
    f3_rows: dict[str, list[dict[str, object]]] = {}
    for kind in ("fold", "corner"):
        f3_rows[kind] = [
            {"id": str(row["id"]), "count": gadget_schedule_overlap(tuple(row["schedule"]), kind, 3)}
            for row in all_grid_rows
        ]
    timings["three_gadget_terms"] = time.monotonic() - t0
    f3_fold_total = sum(int(row["count"]) for row in f3_rows["fold"])
    f3_corner_total = sum(int(row["count"]) for row in f3_rows["corner"])
    three_gadget = gadget_count(R, 3)
    three_gadget_net = three_gadget - f3_fold_total
    expected_f3_nonzero = {
        "H7_9": 31632108085872,
        "H7_11": 51967034712504,
        "H9_7": 31632108085872,
        "H9_9": 424397450152116,
        "H9_11": 166068567450828,
        "H11_7": 51967034712504,
        "H11_9": 166068567450828,
        "H11_11": 9790890598008,
    }
    f3_fold_map = {str(row["id"]): int(row["count"]) for row in f3_rows["fold"]}
    f3_corner_map = {str(row["id"]): int(row["count"]) for row in f3_rows["corner"]}
    check("three_defect_families", three_gadget == 2325336517026900
          and f3_fold_total == f3_corner_total == 933523761248532
          and {key: value for key, value in f3_fold_map.items() if value} == expected_f3_nonzero
          and {key: value for key, value in f3_corner_map.items() if value} == expected_f3_nonzero
          and all(f3_fold_map[identifier] == f3_corner_map[identifier] == 0
                  for identifier, _ in list(VARIABLE_SCHEDULES) + list(THREE_SCHEDULES)),
          "[COMPUTATION] exact F3/Q3 count and all 82 schedule intersections")

    union = old_base + variable_net + three_net + 2 * three_gadget_net
    m_value = C36 - union
    check("union_a35", union == 507841195880595165555,
          "[THEOREM] exact finite inclusion--exclusion gives the enlarged a_35 subfamily")
    check("chain_M", m_value == 2940863015138821131395115,
          "[THEOREM] M=c_36-a_35_lower from the first-backtrack chain")
    check("strict_improvement", union > old_base and m_value < C36 - old_base,
          "[THEOREM] strictly improves the wave-12 lower endpoint without using K_c data")

    lower, upper = atanh_root_interval(m_value)
    floor = floor40(lower)
    check("endpoint_floor", floor == "0.2122129498516328253630821339504383555437",
          "[COMPUTATION] 120-dps directed floor of atanh(M^(-1/36))")
    with localcontext() as context:
        context.prec = 180
        check("endpoint_rounding", Decimal(floor) <= Decimal(lower)
              and Decimal(upper) < Decimal(floor) + Decimal(1).scaleb(-REPORT_PLACES),
              "[COMPUTATION] floor<=interval lower<=true<=upper<floor+10^-40")

    elapsed = time.monotonic() - started
    resource_walls = [
        {
            "kind": "[COMPUTATION] observed",
            "description": "measured producer wall time and peak resident set after all exact DPs",
            "wall_seconds": elapsed,
            "peak_rss_bytes": peak_rss_bytes(),
        },
        {
            "kind": "[UNRESOLVED] preflight",
            "description": "declared input-size wall, not an observed lower/upper resource theorem",
            "profile_max_depth": 13,
            "old_grid_rows": 25,
            "new_grid_rows": 16,
            "pairwise_separator_rows": 820,
            "three_gadget_schedule_intersections": 82,
            "wall_budget_seconds": 300,
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
            "elapsed_seconds": elapsed,
            "timings": timings,
        },
        "data": {
            "classification": "[THEOREM] a_35 >= 507841195880595165555 by an exact finite union; [UNRESOLVED] this is not an all-size SAW-family theorem or a solution of the 3D Ising model.",
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
            },
            "base": {
                "status": "[THEOREM] inherited but independently recomputed",
                "L": l35,
                "C": c35,
                "C_cap_L": c_cap_l,
                "old_height_grid_rows": old_rows,
                "old_grid_net": old_grid_net,
                "F2_count": two_gadget,
                "Q2_count": two_gadget,
                "old_F2_cap_H_sum": old_f2,
                "old_Q2_cap_H_sum": old_q2,
                "old_union": old_base,
            },
            "schedule_disjointness": {
                "status": "[THEOREM]",
                "claim": "All 41 height-slab families are pairwise disjoint; every higher intersection of two or more schedule families is empty because it is contained in a listed empty pair.",
                "method": "for every pair, a listed step has disjoint exact allowed-height sets; each set is bounded by R because no R-step walk can have larger height.",
                "pairwise_separator_rows": separator_rows,
                "higher_intersections": "[THEOREM] all triple and higher schedule intersections are zero by pairwise disjointness; no unlisted nonzero higher schedule term exists.",
            },
            "variable_cut_grid": {
                "status": "[THEOREM]",
                "definition": "V_q=B(13,q), translated B(11,3), translated T(11), q in {5,7,9,11}",
                "rises": [5, 7, 9, 11],
                "injection": "[LEMMA] strict consecutive height bands make each concatenation self-avoiding, fixed cuts recover its components, and positive height avoids e_1.",
                "overlap_method": "[COMPUTATION] signed-orthant Möbius DPs for C∩H; seven regular-language DPs for L∩H and C∩L∩H; two-gadget DPs for F2/Q2.",
                "rows": variable_rows,
                "totals": {"H": variable_h, "C_cap_H": variable_c, "net": variable_net},
            },
            "three_block_grid": {
                "status": "[THEOREM]",
                "definition": "J_(q1,q3)=B(7,q1), translated B(7,3), translated B(9,q3), translated T(12), q1 in {3,5,7}, q3 in {3,5,7,9}",
                "admissible_odd_rises": [3, 5, 7, 9, 11],
                "selected_rise_grid": {"q1": [3, 5, 7], "q2": [3], "q3": [3, 5, 7, 9]},
                "injection": "[LEMMA] four strict height bands are disjoint except at splices, so concatenation is self-avoiding and e_1-avoiding.",
                "overlap_method": "[COMPUTATION] the same exact C/L/F2/Q2 DPs as the variable-cut grid.",
                "rows": three_rows,
                "totals": {"H": three_h, "C_cap_H": three_c, "net": three_net},
            },
            "cm_mobius_tables": {
                "status": "[COMPUTATION]",
                "method": "all signed-coordinate orthant terms for every schedule; sums equal each stored C_cap_H row",
                "rows": mobius_tables,
            },
            "three_defect_families": {
                "fold": {
                    "status": "[THEOREM]",
                    "definition": "u F_y v F_y w F_y x; F_y=(-e1,+e2,+e1), u in M*, v,w in M+ each begin +e2/+e3, x empty or begins +e2/+e3",
                    "count": three_gadget,
                    "self_avoidance": "[LEMMA] contract every F_y to a +e2 skeleton edge; each required post-gadget advance forbids its unique side-vertex return and monotonicity forbids all later collisions.",
                    "e1_avoidance": "[LEMMA] in (a,b,c)=(-x,y,z), the contracted skeleton and each F_y have a>=0, hence x never becomes positive.",
                    "disjointness": "[THEOREM] exactly three +e1 defect letters exclude L and F2; both x signs exclude C; Q2/Q3 contain no +e1.",
                    "cap_all_grid": f3_rows["fold"],
                    "cap_all_grid_sum": f3_fold_total,
                    "net": three_gadget_net,
                },
                "corner": {
                    "status": "[THEOREM]",
                    "definition": "u Q_y v Q_y w Q_y x; Q_y=(+e2,+e3,-e2), u in M*, v,w in M+ each begin -e1/+e3, x empty or begins -e1/+e3",
                    "count": three_gadget,
                    "self_avoidance": "[LEMMA] contract every Q_y to a +e3 skeleton edge; each required post-gadget advance forbids its unique side-vertex return and monotonicity forbids all later collisions.",
                    "e1_avoidance": "[LEMMA] no +e1 letter occurs.",
                    "disjointness": "[THEOREM] exactly three -e2 defect letters exclude L and Q2; both y signs exclude C; F2/F3 contain +e1 whereas Q_y does not.",
                    "cap_all_grid": f3_rows["corner"],
                    "cap_all_grid_sum": f3_corner_total,
                    "net": three_gadget_net,
                },
            },
            "union_certificate": {
                "status": "[THEOREM]",
                "expression": "U_old + sum_newH(H-C∩H-L∩H+C∩L∩H-F2∩H-Q2∩H) + |F3|-sum_H|F3∩H| + |Q3|-sum_H|Q3∩H|",
                "old_union": old_base,
                "variable_cut_net": variable_net,
                "three_block_net": three_net,
                "F3_net": three_gadget_net,
                "Q3_net": three_gadget_net,
                "a35_lower": union,
                "scope": "[UNRESOLVED] finite r=35 certificate only; no all-r lower-family recursion is claimed.",
            },
            "chain": {
                "status": "[THEOREM]",
                "statement": "the previously proved first-backtrack inequality gives mu^n<=c_n-a_(n-1); apply at n=36",
                "M": m_value,
                "mu_power_bound": "mu^36<=2940863015138821131395115",
            },
            "endpoint": {
                "status": "[THEOREM] via the SAW-domination proof in proofs/kc_bounds.md section A",
                "statement": "K_c>=atanh(M^(-1/36))",
                "interval_120dps": [lower, upper],
                "floor_40": floor,
            },
            "resource_walls": resource_walls,
            "unresolved": [
                "[UNRESOLVED] the selected variable-cut and three-block grids are explicit disjoint subfamilies, not an optimization over all 35-step schedules.",
                "[UNRESOLVED] no all-size strengthening of L_r or any Kesten-ratio theorem follows from this finite union.",
                "[UNRESOLVED] the finite certificate does not solve the three-dimensional Ising model.",
            ],
        },
        "checks": checks,
    }
    RESULT.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    print(f"\nwrote {RESULT}")
    print(f"elapsed {elapsed:.1f}s; peak RSS {peak_rss_bytes()} bytes")
    print(f"a_35 lower bound: {union}")
    print(f"M: {m_value}")
    print(f"K_c lower floor (40): {floor}")


if __name__ == "__main__":
    try:
        main()
    except AssertionError as exc:
        print(f"CERTIFICATE FAILURE: {exc}")
        raise
