#!/usr/bin/env python3
"""Clean-room verifier for the e123 variable-cut SAW union certificate.

This verifier deliberately does not import ``experiments/e123_saw_union3.py``.
It rebuilds the positive-height block data by an iterative tuple-coordinate
SAW DFS, recomputes every height-schedule intersection by finite integer DPs,
and checks every one of the 820 cross-grid disjointness certificates.  The
producer artifact is only read after the independent certificate is complete.

Run:
    PYTHONPATH=src .venv/bin/python tests/test_saw_union3.py
"""
from __future__ import annotations

from collections import defaultdict
from decimal import Decimal, ROUND_FLOOR, localcontext
from itertools import product
import hashlib
import json
from pathlib import Path

import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "bounds" / "saw_union3.json"
SAW_SOURCE = ROOT / "sources" / "fulltext" / "schram_barkema_bisseling2011.pdf"

C36 = 2941370856334701726560670
C36_SHA256 = "898d8333e8d70acd279e29f226ff108550b132463aa662f86df590b0ad72cb12"
DPS = 150
R = 35
M = (1, 2, 4)  # -e1, +e2, +e3; every letter raises h=-x+y+z.
STEPS = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))
DH = (-1, 1, 1, -1, 1, -1)

EXPECTED_TAIL = (1, 3, 9, 45, 171, 837, 3411, 16425, 69525, 331947,
                 1436643, 6827355, 29971017, 142016661)
EXPECTED_BLOCKS = {
    1: {1: 3},
    2: {2: 9},
    3: {3: 27},
    4: {4: 81},
    5: {3: 108, 5: 243},
    6: {4: 648, 6: 729},
    7: {3: 432, 5: 2916, 7: 2187},
    8: {4: 5778, 6: 11664, 8: 6561},
    9: {3: 1566, 5: 34668, 7: 43740, 9: 19683},
    10: {4: 48168, 6: 167670, 8: 157464, 10: 59049},
    11: {3: 5832, 5: 411156, 7: 729000, 9: 551124, 11: 177147},
    12: {4: 394524, 6: 2356128, 8: 2969946, 10: 1889568, 12: 531441},
    13: {3: 20952, 5: 4815450, 7: 11544444, 9: 11573604, 11: 6377292,
         13: 1594323},
}

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

EXPECTED_L35 = 644233352324156721
EXPECTED_C35 = 333543290395947705
EXPECTED_CL = 50032713330104219
EXPECTED_OLD_GRID_C = 51027099125602155
EXPECTED_OLD_GRID_L = 491050350045737865
EXPECTED_OLD_GRID_CL = 50032404092458907
EXPECTED_TWO_GADGET = 13710824278006626
EXPECTED_TWO_GADGET_OLD_H = 7732355849776818
EXPECTED_BASE_UNION = 499330428831189067251
EXPECTED_VAR_H = 1366157335427744400
EXPECTED_VAR_C = 4861170126720
EXPECTED_THREE_H = 7141848962632019280
EXPECTED_THREE_C = 18012995095392
EXPECTED_NEW_GRID_NET = 8507983423894541568
EXPECTED_THREE_GADGET = 2325336517026900
EXPECTED_THREE_GADGET_OLD_H = 933523761248532
EXPECTED_THREE_GADGET_NET = 1391812755778368
EXPECTED_UNION = 507841195880595165555
EXPECTED_M = 2940863015138821131395115
EXPECTED_FLOOR = "0.2122129498516328253630821339504383555437"

FAILURES: list[str] = []


def check(name: str, passed: bool, detail: str = "") -> None:
    print(("PASS" if passed else "FAIL") + f": {name}" + (f" ({detail})" if detail else ""))
    if not passed:
        FAILURES.append(name)


# ---------------------------------------------------------------------------
# Independent raw positive-height enumeration.
# ---------------------------------------------------------------------------


def positive_height_profile(max_depth: int) -> tuple[list[int], list[dict[int, int]]]:
    """Iteratively enumerate tuple-coordinate positive-height SAWs.

    A terminal contributes to B(n,q) precisely when its height is a strict
    record: positivity gives the lower interior bound and strictness gives the
    upper bound.  This traversal has no dependency on the producer.
    """
    tail = [0] * (max_depth + 1)
    blocks = [defaultdict(int) for _ in range(max_depth + 1)]
    origin = (0, 0, 0)
    stack: list[tuple[tuple[int, int, int], int, int, int, tuple[tuple[int, int, int], ...]]] = [
        (origin, 0, 0, 0, (origin,))
    ]
    while stack:
        point, height, previous_max, depth, vertices = stack.pop()
        tail[depth] += 1
        if depth and height > previous_max:
            blocks[depth][height] += 1
        if depth == max_depth:
            continue
        x, y, z = point
        next_previous_max = max(previous_max, height)
        for letter, (dx, dy, dz) in enumerate(STEPS):
            nxt = (x + dx, y + dy, z + dz)
            if nxt in vertices:
                continue
            next_height = height + DH[letter]
            if next_height <= 0:
                continue
            stack.append((nxt, next_height, next_previous_max, depth + 1, vertices + (nxt,)))
    return tail, [dict(row) for row in blocks]


# ---------------------------------------------------------------------------
# Generic two-/three-block height schedules and exact word DPs.
# ---------------------------------------------------------------------------


def allowed_heights(schedule: tuple[int, ...], position: int) -> set[int]:
    """Return every allowed height at one-based position, bounded by R.

    The bound is harmless because a 35-step nearest-neighbour word cannot
    attain a height outside [-35,35].  A disjoint returned pair is therefore a
    direct finite certificate that the corresponding walk families are
    disjoint; no heuristic height sampling is involved.
    """
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


def words_in_schedule(letters_by_position: list[tuple[int, ...]], schedule: tuple[int, ...],
                      coordinate_monotone: bool = False) -> int:
    """Exact integer DP for a finite letter language subject to height slabs."""
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


def signed_letters(code: tuple[int | None, int | None, int | None]) -> tuple[int, ...]:
    return tuple(axis * 2 + (0 if sign == 0 else 1)
                 for axis, sign in enumerate(code) if sign is not None)


def cm_overlap_used_patterns(schedule: tuple[int, ...]) -> int:
    """C∩H by used signed-pattern inclusion-exclusion, not producer Möbius."""
    total = 0
    for code in product((None, 0, 1), repeat=3):
        pattern = signed_letters(code)
        if not pattern:
            continue
        exact_pattern = 0
        for mask in range(1 << len(pattern)):
            subset = tuple(pattern[index] for index in range(len(pattern)) if mask & (1 << index))
            exact_pattern += (-1) ** (len(pattern) - len(subset)) * words_in_schedule(
                [subset] * R, schedule
            )
        total += exact_pattern
    return total


# ---------------------------------------------------------------------------
# Existing L and C families, rebuilt without importing their producer.
# ---------------------------------------------------------------------------


def l_count(r: int) -> int:
    def suffix_count(length: int) -> int:
        return 1 if length == 0 else 2 * 3 ** (length - 1)

    corners = sum(3 ** k * suffix_count(r - 2 - k) for k in range(r - 1)) if r >= 2 else 0
    folds = sum(3 ** k * suffix_count(r - 3 - k) for k in range(r - 2)) if r >= 3 else 0
    return 3 ** r + 4 * corners + 2 * folds


def l_schedule_overlap(schedule: tuple[int, ...], coordinate_monotone: bool = False) -> tuple[int, dict[str, int]]:
    """Independently sum every one of the seven disjoint L languages."""
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
                language[defect_position + 1] = tuple(
                    letter for letter in M if letter != forbidden_successor
                )
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


def cm_count(r: int) -> int:
    """Count C by surjective used-sign patterns, independently of its automaton."""
    if r == 0:
        return 1
    total = 0
    for code in product((None, 0, 1), repeat=3):
        pattern = signed_letters(code)
        if not pattern:
            continue
        all_used = 0
        for mask in range(1 << len(pattern)):
            subset_size = sum(1 for index in range(len(pattern)) if mask & (1 << index))
            all_used += (-1) ** (len(pattern) - subset_size) * subset_size ** r
        first_plus = 0
        if 0 in pattern:
            others = tuple(letter for letter in pattern if letter != 0)
            for mask in range(1 << len(others)):
                subset_size = sum(1 for index in range(len(others)) if mask & (1 << index))
                first_plus += (-1) ** (len(others) - subset_size) * (1 + subset_size) ** (r - 1)
        total += all_used - first_plus
    return total


# ---------------------------------------------------------------------------
# Two- and three-gadget regular languages, including direct small SAW checks.
# ---------------------------------------------------------------------------


def gadget_spec(kind: str) -> tuple[tuple[int, int, int], tuple[int, int]]:
    if kind == "fold":
        return (1, 2, 0), (2, 4)       # -e1,+e2,+e1; advance +e2/+e3
    if kind == "corner":
        return (2, 4, 3), (1, 4)       # +e2,+e3,-e2; advance -e1/+e3
    raise ValueError(kind)


def free_segment_lengths(total: int, gadgets: int):
    """Yield n0,n1,...,n_(g-1),ng with internal ni>=1 and final ng>=0."""
    for first in range(total + 1):
        def append_internal(prefix: tuple[int, ...], remaining: int, pending: int):
            if pending == 0:
                yield prefix + (remaining,)
                return
            for length in range(1, remaining - (pending - 1) + 1):
                yield from append_internal(prefix + (length,), remaining - length, pending - 1)
        yield from append_internal((first,), total - first, gadgets - 1)


def protected_segment_words(length: int, safe: tuple[int, int]) -> tuple[tuple[int, ...], ...]:
    if length == 0:
        return ((),)
    return tuple((first,) + rest for first in safe for rest in product(M, repeat=length - 1))


def gadget_words(r: int, kind: str, gadgets: int) -> list[tuple[int, ...]]:
    """Direct finite generator used only for small exact safety checks."""
    total = r - 3 * gadgets
    if total < gadgets - 1:
        return []
    gadget, safe = gadget_spec(kind)
    result: list[tuple[int, ...]] = []
    for lengths in free_segment_lengths(total, gadgets):
        first_length, *rest_lengths = lengths
        internal_lengths = rest_lengths[:-1]
        final_length = rest_lengths[-1]
        choices = [tuple(product(M, repeat=first_length))]
        choices.extend(protected_segment_words(length, safe) for length in internal_lengths)
        choices.append(protected_segment_words(final_length, safe))
        for segments in product(*choices):
            word = segments[0]
            for index in range(gadgets):
                word += gadget
                word += segments[index + 1]
            if len(word) != r:
                raise AssertionError("gadget word length")
            result.append(word)
    return result


def gadget_count(r: int, gadgets: int) -> int:
    total = r - 3 * gadgets
    if total < gadgets - 1:
        return 0
    answer = 0
    for lengths in free_segment_lengths(total, gadgets):
        first_length, *rest_lengths = lengths
        count = 3 ** first_length
        for length in rest_lengths[:-1]:
            count *= 2 * 3 ** (length - 1)
        final_length = rest_lengths[-1]
        count *= 1 if final_length == 0 else 2 * 3 ** (final_length - 1)
        answer += count
    return answer


def gadget_schedule_overlap(schedule: tuple[int, ...], kind: str, gadgets: int) -> int:
    """Exact finite-language DP for one gadget family against one schedule."""
    total = R - 3 * gadgets
    if total < gadgets - 1:
        return 0
    gadget, safe = gadget_spec(kind)
    total_count = 0
    singleton_gadget = [(letter,) for letter in gadget]
    for lengths in free_segment_lengths(total, gadgets):
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
            raise AssertionError("gadget DP language length")
        total_count += words_in_schedule(language, schedule)
    return total_count


def is_e1_avoiding_saw(word: tuple[int, ...]) -> bool:
    point = (0, 0, 0)
    visited = {point}
    for letter in word:
        dx, dy, dz = STEPS[letter]
        point = (point[0] + dx, point[1] + dy, point[2] + dz)
        if point == (1, 0, 0) or point in visited:
            return False
        visited.add(point)
    return True


# ---------------------------------------------------------------------------
# Arithmetic helpers.
# ---------------------------------------------------------------------------


def height_schedule_count(schedule: tuple[int, ...], blocks: list[dict[int, int]], tail: list[int]) -> int:
    if len(schedule) == 5:
        l1, q1, l2, q2, tail_length = schedule
        return blocks[l1][q1] * blocks[l2][q2] * tail[tail_length]
    l1, q1, l2, q2, l3, q3, tail_length = schedule
    return blocks[l1][q1] * blocks[l2][q2] * blocks[l3][q3] * tail[tail_length]


def enclosure(magnitude: int) -> tuple[str, str]:
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
        return format(Decimal(value).quantize(Decimal(1).scaleb(-40), rounding=ROUND_FLOOR), "f")


def main() -> None:
    # RED phase: this fails until the independently verified producer artifact exists.
    if not RESULT.exists():
        raise AssertionError("e123 producer artifact is missing")

    check("external_source_hash", hashlib.sha256(SAW_SOURCE.read_bytes()).hexdigest() == C36_SHA256,
          "cached exact c_36 Table-I source")

    tail, blocks = positive_height_profile(13)
    check("positive_height_tail_profile", tuple(tail) == EXPECTED_TAIL,
          "iterative tuple DFS reproduces T(0)..T(13)")
    check("positive_height_block_profile", blocks == [dict() if index == 0 else EXPECTED_BLOCKS[index]
                                                        for index in range(14)],
          "strict-record terminals reproduce every B(n,q), n<=13")

    l35 = l_count(R)
    c35 = cm_count(R)
    c_cap_l = 3 ** R + 4 * (R - 1) * 2 ** (R - 2)
    check("base_L_C_counts", l35 == EXPECTED_L35 and c35 == EXPECTED_C35 and c_cap_l == EXPECTED_CL,
          "independent seven-language and used-pattern counts")

    # Rebuild all 25 old schedules, including every C/L pair and C/L/H triple.
    old_rows: list[dict[str, object]] = []
    for q1 in OLD_GRID_Q:
        for q2 in OLD_GRID_Q:
            schedule = (11, q1, 11, q2, 13)
            l_overlap, _ = l_schedule_overlap(schedule)
            cl_overlap, _ = l_schedule_overlap(schedule, coordinate_monotone=True)
            old_rows.append({
                "id": f"H{q1}_{q2}",
                "schedule": list(schedule),
                "H": height_schedule_count(schedule, blocks, tail),
                "C_cap_H": cm_overlap_used_patterns(schedule),
                "L_cap_H": l_overlap,
                "C_cap_L_cap_H": cl_overlap,
                "F2_cap_H": gadget_schedule_overlap(schedule, "fold", 2),
                "Q2_cap_H": gadget_schedule_overlap(schedule, "corner", 2),
            })
    old_grid_net = sum(int(row["H"]) - int(row["C_cap_H"]) - int(row["L_cap_H"])
                       + int(row["C_cap_L_cap_H"]) for row in old_rows)
    check("old_grid_all_overlap_terms", len(old_rows) == 25
          and sum(int(row["C_cap_H"]) for row in old_rows) == EXPECTED_OLD_GRID_C
          and sum(int(row["L_cap_H"]) for row in old_rows) == EXPECTED_OLD_GRID_L
          and sum(int(row["C_cap_L_cap_H"]) for row in old_rows) == EXPECTED_OLD_GRID_CL
          and sum(int(row["F2_cap_H"]) for row in old_rows) == EXPECTED_TWO_GADGET_OLD_H
          and sum(int(row["Q2_cap_H"]) for row in old_rows) == EXPECTED_TWO_GADGET_OLD_H,
          "all 25 old C/L/F2/Q2 pair and triple schedule intersections")

    # Directly generate small languages before trusting their general count formula.
    small_languages_ok = True
    for gadgets, lengths in ((2, range(8, 12)), (3, range(11, 15))):
        for length in lengths:
            fold_words = gadget_words(length, "fold", gadgets)
            corner_words = gadget_words(length, "corner", gadgets)
            if (len(fold_words) != gadget_count(length, gadgets)
                    or len(corner_words) != gadget_count(length, gadgets)
                    or not all(is_e1_avoiding_saw(word) for word in fold_words + corner_words)
                    or set(fold_words).intersection(corner_words)):
                small_languages_ok = False
    check("two_and_three_gadget_small_saw_certificate", small_languages_ok,
          "direct words catch unsafe successor, sign, and gadget-count defects")

    two_gadget = gadget_count(R, 2)
    old_base = l35 + c35 - c_cap_l + old_grid_net + 2 * (two_gadget - EXPECTED_TWO_GADGET_OLD_H)
    check("old_union_rebuilt", two_gadget == EXPECTED_TWO_GADGET and old_base == EXPECTED_BASE_UNION,
          "the wave-12 finite union is recomputed from raw profile and DPs")

    # Every selected new schedule has a direct height separator from every old
    # or new schedule.  This lists ALL pairwise H_i∩H_j=empty facts; every
    # triple/higher grid intersection is then a subset of a recorded zero pair.
    new_specs = list(VARIABLE_SCHEDULES) + list(THREE_SCHEDULES)
    schedule_specs = [(str(row["id"]), tuple(row["schedule"])) for row in old_rows] + list(new_specs)
    separator_rows: list[dict[str, object]] = []
    for left_index, (left_id, left_schedule) in enumerate(schedule_specs):
        for right_id, right_schedule in schedule_specs[left_index + 1:]:
            separator = first_separator(left_schedule, right_schedule)
            if separator is None:
                raise AssertionError(f"no direct separator for {left_id},{right_id}")
            separator_rows.append({"left": left_id, "right": right_id, **separator})
    check("all_grid_pairs_directly_separated", len(schedule_specs) == 41
          and len(separator_rows) == 41 * 40 // 2,
          "each of the 820 pairs has an exact disjoint allowed-height step")

    def build_new_rows(category: str, specs: tuple[tuple[str, tuple[int, ...]], ...]) -> list[dict[str, object]]:
        rows: list[dict[str, object]] = []
        for identifier, schedule in specs:
            l_overlap, _ = l_schedule_overlap(schedule)
            cl_overlap, _ = l_schedule_overlap(schedule, coordinate_monotone=True)
            rows.append({
                "id": identifier,
                "schedule": list(schedule),
                "H": height_schedule_count(schedule, blocks, tail),
                "C_cap_H": cm_overlap_used_patterns(schedule),
                "L_cap_H": l_overlap,
                "C_cap_L_cap_H": cl_overlap,
                "F2_cap_H": gadget_schedule_overlap(schedule, "fold", 2),
                "Q2_cap_H": gadget_schedule_overlap(schedule, "corner", 2),
                "category": category,
            })
        return rows

    variable_rows = build_new_rows("variable_cut", VARIABLE_SCHEDULES)
    three_rows = build_new_rows("three_block", THREE_SCHEDULES)
    new_rows = variable_rows + three_rows
    variable_net = sum(int(row["H"]) - int(row["C_cap_H"]) - int(row["L_cap_H"])
                       + int(row["C_cap_L_cap_H"]) - int(row["F2_cap_H"]) - int(row["Q2_cap_H"])
                       for row in variable_rows)
    three_net = sum(int(row["H"]) - int(row["C_cap_H"]) - int(row["L_cap_H"])
                    + int(row["C_cap_L_cap_H"]) - int(row["F2_cap_H"]) - int(row["Q2_cap_H"])
                    for row in three_rows)
    check("variable_cut_grid_every_term", len(variable_rows) == 4
          and sum(int(row["H"]) for row in variable_rows) == EXPECTED_VAR_H
          and sum(int(row["C_cap_H"]) for row in variable_rows) == EXPECTED_VAR_C
          and all(int(row[key]) == 0 for row in variable_rows
                  for key in ("L_cap_H", "C_cap_L_cap_H", "F2_cap_H", "Q2_cap_H")),
          "all four variable-cut pair/triple terms are exact, including zeros")
    check("three_block_grid_every_term", len(three_rows) == 12
          and sum(int(row["H"]) for row in three_rows) == EXPECTED_THREE_H
          and sum(int(row["C_cap_H"]) for row in three_rows) == EXPECTED_THREE_C
          and all(int(row[key]) == 0 for row in three_rows
                  for key in ("L_cap_H", "C_cap_L_cap_H", "F2_cap_H", "Q2_cap_H")),
          "all twelve three-block pair/triple terms are exact, including zeros")
    check("new_grid_net", variable_net + three_net == EXPECTED_NEW_GRID_NET,
          "exact inclusion-exclusion relative to the full old union")

    # The two new three-gadget languages have exact intersections with all 41
    # disjoint height schedules.  Their only nonzero rows occur in the old grid.
    all_grid_rows = old_rows + new_rows
    f3_rows: dict[str, list[dict[str, object]]] = {}
    for kind in ("fold", "corner"):
        f3_rows[kind] = [
            {"id": str(row["id"]), "count": gadget_schedule_overlap(tuple(row["schedule"]), kind, 3)}
            for row in all_grid_rows
        ]
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
    three_gadget = gadget_count(R, 3)
    fold_f3 = {str(row["id"]): int(row["count"]) for row in f3_rows["fold"]}
    corner_f3 = {str(row["id"]): int(row["count"]) for row in f3_rows["corner"]}
    check("three_gadget_every_grid_intersection", three_gadget == EXPECTED_THREE_GADGET
          and sum(fold_f3.values()) == sum(corner_f3.values()) == EXPECTED_THREE_GADGET_OLD_H
          and {key: value for key, value in fold_f3.items() if value} == expected_f3_nonzero
          and {key: value for key, value in corner_f3.items() if value} == expected_f3_nonzero
          and all(fold_f3[identifier] == corner_f3[identifier] == 0
                  for identifier, _ in new_specs),
          "all 41 F3/Q3 schedule terms; 16 new-grid zero terms retained")

    three_gadget_net = three_gadget - sum(fold_f3.values())
    union = old_base + variable_net + three_net + 2 * three_gadget_net
    m_value = C36 - union
    check("full_inclusion_exclusion", three_gadget_net == EXPECTED_THREE_GADGET_NET
          and union == EXPECTED_UNION and m_value == EXPECTED_M,
          "base plus variable/three-block grids plus F3/Q3 families")
    no_new_c_overlap = old_base + sum(int(row["H"]) for row in new_rows) + 2 * three_gadget_net
    no_f3_overlap = old_base + variable_net + three_net + 2 * three_gadget
    check("overlap_omission_traps", no_new_c_overlap - union == EXPECTED_VAR_C + EXPECTED_THREE_C
          and no_f3_overlap - union == 2 * EXPECTED_THREE_GADGET_OLD_H,
          "dropping a new C∩H or F3∩H term changes the asserted union")
    check("strict_improvement", union > EXPECTED_BASE_UNION and m_value < C36 - EXPECTED_BASE_UNION,
          "strictly improves the prior finite a_35 certificate")

    lo, hi = enclosure(m_value)
    floor = floor40(lo)
    check("endpoint_floor_150dps", floor == EXPECTED_FLOOR,
          "independent 150-dps directed lower rounding")
    with localcontext() as context:
        context.prec = 220
        check("endpoint_directed_semantics", Decimal(floor) <= Decimal(lo)
              and Decimal(hi) < Decimal(floor) + Decimal(1).scaleb(-40),
              "floor <= interval lower <= true endpoint <= upper < floor+10^-40")

    payload = json.loads(RESULT.read_text(encoding="utf-8"))
    data = payload["data"]
    check("artifact_provenance", payload["provenance"]["script"] == "experiments/e123_saw_union3.py"
          and payload["provenance"]["mpmath_iv_dps"] >= 120
          and all(row["passed"] for row in payload["checks"]),
          "producer provenance and checks")
    check("artifact_base_rows", data["base"]["old_height_grid_rows"] == old_rows,
          "every inherited C/L/F2/Q2 grid term equals clean-room DP")
    check("artifact_disjointness_rows", data["schedule_disjointness"]["pairwise_separator_rows"] == separator_rows,
          "all 820 direct height separators are retained verbatim")
    check("artifact_new_grid_rows", data["variable_cut_grid"]["rows"] == variable_rows
          and data["three_block_grid"]["rows"] == three_rows,
          "every new pair and triple schedule term equals clean-room DP")
    check("artifact_three_gadget_rows", data["three_defect_families"]["fold"]["cap_all_grid"] == f3_rows["fold"]
          and data["three_defect_families"]["corner"]["cap_all_grid"] == f3_rows["corner"],
          "all 82 F3/Q3 schedule intersections equal clean-room DP")
    check("artifact_final_values", data["union_certificate"]["a35_lower"] == union
          and data["chain"]["M"] == m_value
          and data["endpoint"]["floor_40"] == floor
          and any(row.get("kind") == "[COMPUTATION] observed" for row in data["resource_walls"])
          and any(row.get("kind") == "[UNRESOLVED] preflight" for row in data["resource_walls"]),
          "stored union, endpoint, and observed/preflight resource distinctions")

    print()
    if FAILURES:
        print(f"FAILED ({len(FAILURES)}): {FAILURES}")
        raise SystemExit(1)
    print("OK")


if __name__ == "__main__":
    main()
