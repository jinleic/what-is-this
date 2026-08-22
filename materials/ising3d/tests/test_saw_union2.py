#!/usr/bin/env python3
"""Clean-room verifier for the e120 SAW-union extension.

This file intentionally does not import ``experiments/e120_saw_union2.py``.
It re-enumerates the positive-height profile by an iterative tuple-coordinate
DFS, uses used-sign-pattern inclusion--exclusion for coordinate-monotone
intersections, and uses independent finite word DPs for the regular defect
languages.  It also contains overlap-omission traps and a 150-dps directed
rounding check.

Run: PYTHONPATH=src .venv/bin/python tests/test_saw_union2.py
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
RESULT = ROOT / "results" / "bounds" / "saw_union2.json"
SAW_SOURCE = ROOT / "sources" / "fulltext" / "schram_barkema_bisseling2011.pdf"

C36 = 2941370856334701726560670
C36_SHA256 = "898d8333e8d70acd279e29f226ff108550b132463aa662f86df590b0ad72cb12"
DPS = 150
R = 35
M = (1, 2, 4)  # -e1, +e2, +e3
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

EXPECTED_L35 = 644233352324156721
EXPECTED_C35 = 333543290395947705
EXPECTED_CL = 50032713330104219
EXPECTED_GRID_C_SUM = 51027099125602155
EXPECTED_GRID_L_SUM = 491050350045737865
EXPECTED_GRID_CL_SUM = 50032404092458907
EXPECTED_GADGET = 13710824278006626
EXPECTED_GADGET_H_SUM = 7732355849776818
EXPECTED_UNION = 499330428831189067251
EXPECTED_M = 2940871525905870537493419
EXPECTED_FLOOR = "0.2122129322754723621039731646196408804906"

FAILURES: list[str] = []


def check(name: str, passed: bool, detail: str = "") -> None:
    print(("PASS" if passed else "FAIL") + f": {name}" + (f" ({detail})" if detail else ""))
    if not passed:
        FAILURES.append(name)


# ---------------------------------------------------------------------------
# Independent tuple-coordinate profile enumeration.
# ---------------------------------------------------------------------------


def positive_height_profile(max_depth: int) -> tuple[list[int], list[dict[int, int]]]:
    """Enumerate all positive-height SAWs through ``max_depth`` iteratively.

    A terminal at height q contributes to B(n,q) exactly when q is a strict
    record height: then every preceding non-origin height lies in (0,q).
    This reconstructs every B(n,q) and T(n) in one tuple-coordinate traversal.
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
# Height-word DPs for coordinate-monotone and regular-language intersections.
# ---------------------------------------------------------------------------


def schedule_allows(height: int, position: int, schedule: tuple[int, int, int, int, int]) -> bool:
    """Check the positive-band condition after the one-based step ``position``."""
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


def words_in_schedule(letters_by_position: list[tuple[int, ...]], schedule: tuple[int, int, int, int, int],
                      coordinate_monotone: bool = False) -> int:
    """Exact height DP for a prescribed finite step language.

    When ``coordinate_monotone`` is true, the sign slots are carried in the
    state.  No visited set is needed in that case, because a signed-coordinate
    word is self-avoiding.
    """
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

    dp_height: dict[int, int] = {0: 1}
    for position, letters in enumerate(letters_by_position, start=1):
        next_dp_height: dict[int, int] = defaultdict(int)
        for height, count in dp_height.items():
            for letter in letters:
                next_height = height + DH[letter]
                if schedule_allows(next_height, position, schedule):
                    next_dp_height[next_height] += count
        dp_height = next_dp_height
    return sum(dp_height.values())


def signed_letters(code: tuple[int | None, int | None, int | None]) -> tuple[int, ...]:
    return tuple(axis * 2 + (0 if sign == 0 else 1)
                 for axis, sign in enumerate(code) if sign is not None)


def cm_overlap_used_patterns(schedule: tuple[int, int, int, int, int]) -> int:
    """C cap H by exact used-sign-pattern inclusion--exclusion.

    This is deliberately not the producer's global orthant Möbius aggregation:
    it first counts the words using *exactly* each signed set P by inversion
    over the subsets Q of P, then sums over P.
    """
    total = 0
    for code in product((None, 0, 1), repeat=3):
        pattern = signed_letters(code)
        if not pattern:
            continue
        exact_pattern = 0
        for mask in range(1 << len(pattern)):
            subset = tuple(pattern[index] for index in range(len(pattern)) if mask & (1 << index))
            value = words_in_schedule([subset] * R, schedule, coordinate_monotone=False)
            exact_pattern += (-1) ** (len(pattern) - len(subset)) * value
        total += exact_pattern
    return total


# ---------------------------------------------------------------------------
# The L-family and its intersections, independently expressed as fixed-word
# language DPs rather than producer automata.
# ---------------------------------------------------------------------------


def l_count(r: int) -> int:
    def suffix_count(length: int) -> int:
        return 1 if length == 0 else 2 * 3 ** (length - 1)

    corners = sum(3 ** k * suffix_count(r - 2 - k) for k in range(r - 1)) if r >= 2 else 0
    folds = sum(3 ** k * suffix_count(r - 3 - k) for k in range(r - 2)) if r >= 3 else 0
    return 3 ** r + 4 * corners + 2 * folds


def l_schedule_overlap(schedule: tuple[int, int, int, int, int], coordinate_monotone: bool = False) -> tuple[int, dict[str, int]]:
    """Sum the seven disjoint L languages using explicit defect locations."""
    total_length = sum((schedule[0], schedule[2], schedule[4]))
    assert total_length == R
    rows: dict[str, int] = {}
    rows["mono"] = words_in_schedule([M] * R, schedule, coordinate_monotone)

    # (name, defect, required predecessor, forbidden successor)
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
# Candidate-language construction and direct small-length safety checks.
# ---------------------------------------------------------------------------


def gadget_words(r: int, kind: str) -> list[tuple[int, ...]]:
    """Generate the two-gadget language directly for small r.

    ``fold`` uses (-e1,+e2,+e1), while ``corner`` uses
    (+e2,+e3,-e2).  In both cases the inter-gadget and optional final
    suffix must begin with a safe advance, preventing a local revisit.
    """
    if r < 7:
        return []
    if kind == "fold":
        gadget = (1, 2, 0)
        safe = (2, 4)
    elif kind == "corner":
        gadget = (2, 4, 3)
        safe = (1, 4)
    else:
        raise ValueError(kind)
    result: list[tuple[int, ...]] = []
    free_total = r - 6
    for first_length in range(free_total + 1):
        for middle_length in range(1, free_total - first_length + 1):
            last_length = free_total - first_length - middle_length
            for first in product(M, repeat=first_length):
                for middle_rest in product(M, repeat=middle_length - 1):
                    for middle_start in safe:
                        middle = (middle_start,) + middle_rest
                        if last_length == 0:
                            result.append(first + gadget + middle + gadget)
                        else:
                            for last_rest in product(M, repeat=last_length - 1):
                                for last_start in safe:
                                    result.append(first + gadget + middle + gadget + (last_start,) + last_rest)
    return result


def gadget_count(r: int) -> int:
    """Exact sum over the three free-segment lengths, valid for every r>=7."""
    if r < 7:
        return 0
    free_total = r - 6
    total = 0
    for first_length in range(free_total + 1):
        for middle_length in range(1, free_total - first_length + 1):
            last_length = free_total - first_length - middle_length
            total += 3 ** first_length * (2 * 3 ** (middle_length - 1))
            if last_length:
                total += 0  # the preceding line is replaced below for the nonempty tail case
                total -= 3 ** first_length * (2 * 3 ** (middle_length - 1))
                total += 3 ** first_length * (2 * 3 ** (middle_length - 1)) * (2 * 3 ** (last_length - 1))
    return total


def gadget_formula_35() -> int:
    # For N=R-6=29: 2N*3^(N-1) + 2N(N-1)*3^(N-2).
    n = R - 6
    return 2 * n * 3 ** (n - 1) + 2 * n * (n - 1) * 3 ** (n - 2)


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


def gadget_schedule_overlap(schedule: tuple[int, int, int, int, int], kind: str) -> int:
    if kind == "fold":
        gadget = ((1,), (2,), (0,))
        safe = (2, 4)
    elif kind == "corner":
        gadget = ((2,), (4,), (3,))
        safe = (1, 4)
    else:
        raise ValueError(kind)
    total = 0
    free_total = R - 6
    for first_length in range(free_total + 1):
        for middle_length in range(1, free_total - first_length + 1):
            last_length = free_total - first_length - middle_length
            language: list[tuple[int, ...]] = [M] * first_length
            language.extend(gadget)
            language.append(safe)
            language.extend([M] * (middle_length - 1))
            language.extend(gadget)
            if last_length:
                language.append(safe)
                language.extend([M] * (last_length - 1))
            assert len(language) == R
            total += words_in_schedule(language, schedule, coordinate_monotone=False)
    return total


def one_defect_schedule_overlap(schedule: tuple[int, int, int, int, int], defect: int, free: tuple[int, int]) -> int:
    total = 0
    for position in range(R):
        language = [free] * R
        language[position] = (defect,)
        total += words_in_schedule(language, schedule, coordinate_monotone=True)
    return total


# ---------------------------------------------------------------------------
# Coordinate-monotone count, by used-pattern surjection inclusion-exclusion.
# ---------------------------------------------------------------------------


def cm_count(r: int) -> int:
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
# Search tables and interval helpers.
# ---------------------------------------------------------------------------


def two_block_table(blocks: list[dict[int, int]], tail: list[int]) -> list[dict[str, int]]:
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
                            rows.append({"l1": l1, "q1": q1, "l2": l2, "q2": q2,
                                         "l3": l3, "q3": q3, "t": tail_length,
                                         "B1": b1, "B2": b2, "B3": b3, "T": tail[tail_length],
                                         "product": b1 * b2 * b3 * tail[tail_length]})
    return rows


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
    # This deliberately fails before the expensive independent enumeration in
    # the RED phase, and keeps the verifier independent of the producer code.
    if not RESULT.exists():
        raise AssertionError("e120 producer artifact is missing")

    check("external_source_hash", hashlib.sha256(SAW_SOURCE.read_bytes()).hexdigest() == C36_SHA256,
          "cached c_36 source hash")

    tail, blocks = positive_height_profile(13)
    check("positive_height_tail_profile", tuple(tail) == EXPECTED_TAIL,
          "tuple-coordinate iterative DFS reproduces T(0)..T(13)")
    check("positive_height_block_profile", blocks == [dict() if index == 0 else EXPECTED_BLOCKS[index]
                                                        for index in range(14)],
          "strict-record terminals reproduce every B(n,q), n<=13")

    two_rows = two_block_table(blocks, tail)
    three_rows = three_block_table(blocks, tail)
    two_max = max(row["product"] for row in two_rows)
    three_max = max(row["product"] for row in three_rows)
    check("two_block_full_search", len(two_rows) == 392
          and two_max == 75473476338501000000
          and sum(row["product"] == two_max for row in two_rows) == 1,
          "392 exact schedules; unique maximum is (11,7,11,7,13)")
    check("three_block_full_search", len(three_rows) == 1651
          and three_max == 17606452560245513280
          and sum(row["product"] == three_max for row in three_rows) == 6,
          "1651 exact schedules; six permutation-equivalent maxima")

    l35 = l_count(R)
    c35 = cm_count(R)
    c_cap_l = 3 ** R + 4 * (R - 1) * 2 ** (R - 2)
    check("L35_closed_form", l35 == EXPECTED_L35, "seven-family defect lower family")
    check("C35_surjective_pattern_count", c35 == EXPECTED_C35,
          "coordinate-monotone count by per-pattern surjections")
    check("C_cap_L_formula", c_cap_l == EXPECTED_CL, "mono plus four corner-overlap languages")

    grid_rows: list[dict[str, int]] = []
    q_values = sorted(blocks[11])
    for q1 in q_values:
        for q2 in q_values:
            schedule = (11, q1, 11, q2, 13)
            h_count = blocks[11][q1] * blocks[11][q2] * tail[13]
            c_overlap = cm_overlap_used_patterns(schedule)
            l_overlap, l_parts = l_schedule_overlap(schedule, coordinate_monotone=False)
            cl_overlap, cl_parts = l_schedule_overlap(schedule, coordinate_monotone=True)
            grid_rows.append({"q1": q1, "q2": q2, "H": h_count, "C_cap_H": c_overlap,
                              "L_cap_H": l_overlap, "C_cap_L_cap_H": cl_overlap,
                              "L_parts": l_parts, "CL_parts": cl_parts})
    check("height_grid_has_25_disjoint_schedules", len(grid_rows) == 25,
          "different (h_11,h_22) pairs make the height schedules pairwise disjoint")
    check("grid_overlap_sums", sum(row["C_cap_H"] for row in grid_rows) == EXPECTED_GRID_C_SUM
          and sum(row["L_cap_H"] for row in grid_rows) == EXPECTED_GRID_L_SUM
          and sum(row["C_cap_L_cap_H"] for row in grid_rows) == EXPECTED_GRID_CL_SUM,
          "all C/L pairwise and triple intersections recomputed")
    nonzero_l = {(row["q1"], row["q2"]): row["L_cap_H"] for row in grid_rows if row["L_cap_H"]}
    nonzero_cl = {(row["q1"], row["q2"]): row["C_cap_L_cap_H"] for row in grid_rows if row["C_cap_L_cap_H"]}
    check("grid_high_deficit_intersections", nonzero_l == {
              (9, 11): 121063985671653612,
              (11, 9): 121063985671653612,
              (11, 11): 248922378702430641,
          } and nonzero_cl == {
              (9, 11): 240518168576,
              (11, 9): 240518168576,
              (11, 11): 50031923056121755,
          }, "only deficit 0 or 2 schedules can meet a one-defect L word")

    # The prompt's broad one-defect candidate classes are already subsets of C.
    r_one_defect = R * 2 ** (R - 1)
    ry_h = { (row["q1"], row["q2"]): one_defect_schedule_overlap(
        (11, row["q1"], 11, row["q2"], 13), 3, (1, 4)) for row in grid_rows }
    rz_h = { (row["q1"], row["q2"]): one_defect_schedule_overlap(
        (11, row["q1"], 11, row["q2"], 13), 5, (1, 2)) for row in grid_rows }
    expected_one_h = {(9, 11): 120259084288, (11, 9): 120259084288, (11, 11): 188978561024}
    check("one_defect_classes_are_C_redundant", r_one_defect == 35 * 2 ** 34
          and sum(ry_h.values()) == sum(rz_h.values()) == 429496729600
          and {key: value for key, value in ry_h.items() if value} == expected_one_h
          and {key: value for key, value in rz_h.items() if value} == expected_one_h,
          "R_y,R_z subset C; their exact H-grid intersections are retained for audit")

    # Two local gadgets: a two-fold and a two-corner class.  Small direct word
    # generation tests the actual SAW/e1 claim, independently of their formulas.
    small_gadgets_ok = True
    for r in range(7, 11):
        fold_words = gadget_words(r, "fold")
        corner_words = gadget_words(r, "corner")
        if (len(fold_words) != gadget_count(r) or len(corner_words) != gadget_count(r)
                or not all(is_e1_avoiding_saw(word) for word in fold_words + corner_words)
                or set(fold_words).intersection(corner_words)):
            small_gadgets_ok = False
    check("two_gadget_small_saw_certificate", small_gadgets_ok,
          "direct tuples prove both prescribed local languages self-avoiding/e1-avoiding through r=10")
    gadget35 = gadget_formula_35()
    fold_h = {(row["q1"], row["q2"]): gadget_schedule_overlap(
        (11, row["q1"], 11, row["q2"], 13), "fold") for row in grid_rows}
    corner_h = {(row["q1"], row["q2"]): gadget_schedule_overlap(
        (11, row["q1"], 11, row["q2"], 13), "corner") for row in grid_rows}
    expected_gadget_h = {
        (7, 11): 183014339639688,
        (9, 9): 1494617107057452,
        (9, 11): 2455442390165814,
        (11, 7): 183014339639688,
        (11, 9): 2455442390165814,
        (11, 11): 960825283108362,
    }
    check("two_gadget_exact_counts", gadget35 == EXPECTED_GADGET
          and sum(fold_h.values()) == sum(corner_h.values()) == EXPECTED_GADGET_H_SUM
          and {key: value for key, value in fold_h.items() if value} == expected_gadget_h
          and {key: value for key, value in corner_h.items() if value} == expected_gadget_h,
          "both two-defect languages have exact grid intersections")

    # A single safe y-fold is exactly the existing D1 language, not a new family.
    d1_h = {(row["q1"], row["q2"]): row["L_parts"]["D1"] for row in grid_rows}
    check("one_fold_is_L_redundant", sum(d1_h.values()) == 31501343210481297
          and {key: value for key, value in d1_h.items() if value} == {
              (9, 11): 8647427547975258,
              (11, 9): 8647427547975258,
              (11, 11): 14206488114530781,
          }, "the explicitly audited one-fold language is D1 subset L, so it adds no union mass")

    grid_net = sum(row["H"] - row["C_cap_H"] - row["L_cap_H"] + row["C_cap_L_cap_H"]
                   for row in grid_rows)
    fold_net = gadget35 - sum(fold_h.values())
    corner_net = gadget35 - sum(corner_h.values())
    union = l35 + c35 - c_cap_l + grid_net + fold_net + corner_net
    m_value = C36 - union
    check("inclusion_exclusion_identity", union == EXPECTED_UNION and m_value == EXPECTED_M,
          "exact five-part union and c_36 subtraction")
    naive_without_intersections = l35 + c35 + sum(row["H"] for row in grid_rows) + 2 * gadget35
    no_triple = l35 + c35 - c_cap_l + sum(
        row["H"] - row["C_cap_H"] - row["L_cap_H"] for row in grid_rows) + fold_net + corner_net
    no_gadget_overlap = l35 + c35 - c_cap_l + grid_net + 2 * gadget35
    check("trap_every_overlap_class_matters", naive_without_intersections != union
          and no_triple != union and union - no_triple == EXPECTED_GRID_CL_SUM
          and no_gadget_overlap != union and no_gadget_overlap - union == 2 * EXPECTED_GADGET_H_SUM,
          "omitting pair or triple overlap terms changes the claimed lower family")
    check("strict_improvement", union > 76401062626946721319 and m_value < 2941294455272074779839351,
          "strictly stronger than the wave-11 H union L union C certificate")

    payload = json.loads(RESULT.read_text(encoding="utf-8"))
    data = payload["data"]
    check("json_provenance", payload["provenance"]["script"] == "experiments/e120_saw_union2.py"
          and payload["provenance"]["mpmath_iv_dps"] >= 120
          and all(row["passed"] for row in payload["checks"]), "artifact provenance and producer checks")
    check("json_full_search_tables", data["split_search"]["two_block_table"] == two_rows
          and data["split_search"]["three_block_table"] == three_rows,
          "stored tables exactly equal independently recomputed schedules")
    json_grid = data["height_grid"]["rows"]
    check("json_grid_constants", len(json_grid) == len(grid_rows)
          and all(json_row["H"] == independent["H"]
                  and json_row["C_cap_H"] == independent["C_cap_H"]
                  and json_row["L_cap_H"] == independent["L_cap_H"]
                  and json_row["C_cap_L_cap_H"] == independent["C_cap_L_cap_H"]
                  for json_row, independent in zip(json_grid, grid_rows)),
          "stored all-pair and triple grid intersections equal independent values")
    check("json_final_constants", data["union_certificate"]["a35_lower"] == union
          and data["chain"]["M"] == m_value
          and data["two_defect_families"]["fold"]["count"] == gadget35
          and data["two_defect_families"]["corner"]["count"] == gadget35,
          "stored union, base, and two-defect counts")

    lo, hi = enclosure(m_value)
    floor = floor40(lo)
    check("endpoint_floor_150dps", floor == EXPECTED_FLOOR,
          "directed 150-dps downward floor")
    with localcontext() as context:
        context.prec = 220
        check("endpoint_directed_semantics", Decimal(floor) <= Decimal(lo)
              and Decimal(hi) < Decimal(floor) + Decimal(1).scaleb(-40),
              "floor <= interval lower <= true value <= interval upper < floor+10^-40")
    check("json_endpoint", data["endpoint"]["floor_40"] == floor
          and data["endpoint"]["interval_120dps"][0] <= lo <= data["endpoint"]["interval_120dps"][1],
          "stored producer enclosure contains the independent lower endpoint")

    print()
    if FAILURES:
        print(f"FAILED ({len(FAILURES)}): {FAILURES}")
        raise SystemExit(1)
    print("OK")


if __name__ == "__main__":
    main()
