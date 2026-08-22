#!/usr/bin/env python3
"""Clean-room standalone verifier for ``experiments/e92_saw_refine``.

This verifier never imports the producer.  It re-enumerates every small
walk class from scratch (c, a, b, turn, fold, avoid-turn series to depth 10),
checks finite instances of the injection identity exhaustively through m+n=9,
and recomputes every claimed integer certificate -- including the strict
improvement of the K_c lower endpoint -- from hash-audited external constants.
The all-size statements are established by the geometric proofs in the proof
note and by the Fekete limit argument, not by these finite checks.
"""
from __future__ import annotations

from decimal import Decimal, ROUND_FLOOR, localcontext
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]

import mpmath as mp  # noqa: E402


RESULT = ROOT / "results" / "bounds" / "saw_refine.json"
SAW_SOURCE = ROOT / "sources" / "fulltext" / "schram_barkema_bisseling2011.pdf"
C36 = 2941370856334701726560670
C36_SHA256 = "898d8333e8d70acd279e29f226ff108550b132463aa662f86df590b0ad72cb12"
INCUMBENT_K_LOWER = "0.2122119011661678393310862783954278184914"
EXTERNAL_COUNTS = {
    30: 270569905525454674614,
    31: 1274191064726416905966,
    32: 5997359460809616886494,
    33: 28233744272563685150118,
    34: 132853629626823234210582,
    35: 625248129452557974777990,
    36: C36,
}
STORED_C = (1, 6, 30, 150, 726, 3534, 16926, 81390, 387966, 1853886, 8809878)
STORED_A = (1, 5, 25, 121, 589, 2821, 13565, 64661, 308981, 1468313, 6989025)
DPS = 90
DEPTH = 10
PAIR_TOTAL_MAX = 9
FAILURES: list[str] = []

# Expected series re-enumerated independently by this file's own DFS.
EXPECTED_B = (1, 5, 24, 117, 559, 2690, 12817, 61263, 291093, 1385747, 6570755)
EXPECTED_T = (0, 0, 24, 120, 576, 2808, 13416, 64560, 307608, 1470312, 6986232)
EXPECTED_D = (0, 0, 0, 24, 96, 480, 2256, 10944, 51888, 249024, 1181184)
EXPECTED_ATURN = (0, 0, 20, 96, 468, 2236, 10760, 51268, 245052, 1164372, 5542988)
HEIGHT_BLOCK_11_7 = 729000
HALFSPACE_13 = 142016661
HEIGHT_LOWER_35 = 75473476338501000000

OFF = 32
STEPS = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))
AXIS = (0, 0, 1, 1, 2, 2)
SIGN = (1, -1, 1, -1, 1, -1)
MONO_STEPS = (1, 2, 4)


def pack(x: int, y: int, z: int) -> int:
    return ((x + OFF) << 14) | ((y + OFF) << 7) | (z + OFF)


def pack_delta(dx: int, dy: int, dz: int) -> int:
    return (dx << 14) + (dy << 7) + dz


PACKED_STEPS = tuple(pack_delta(*d) for d in STEPS)
FORBID_A = pack(1, 0, 0)
FORBID_B = (pack(1, 0, 0), pack(1, 1, 0))


def check(name: str, passed: bool, detail: str = "") -> None:
    print(("PASS" if passed else "FAIL") + f": {name}" + (f" ({detail})" if detail else ""))
    if not passed:
        FAILURES.append(name)


def enumerate_walks(max_depth: int, forbidden: tuple[int, ...] = (), classify: bool = False):
    result = {key: [0] * (max_depth + 1) for key in ("counts", "turn", "fold")}
    for name in ("mono", "G1", "G2", "G3", "G4", "D1", "D2"):
        result[name] = [0] * (max_depth + 1)
    result["mono"][0] = 1
    forbid = frozenset(forbidden)
    start = pack(0, 0, 0)
    steps_hist: list[int] = []

    def classify_string(depth: int) -> None:
        steps = steps_hist
        for s in steps:
            if s not in MONO_STEPS:
                break
        else:
            result["mono"][depth] += 1
            return
        t = next(i for i, s in enumerate(steps) if s not in MONO_STEPS)
        if any(s not in MONO_STEPS for s in steps[t + 1:]):
            return
        bad = steps[t]
        prev = steps[t - 1] if t >= 1 else None
        nxt = steps[t + 1] if t + 1 < depth else None
        if bad == 3:
            if prev == 4 and nxt != 2:
                result["G1"][depth] += 1
            elif prev == 1 and nxt != 2:
                result["G3"][depth] += 1
        elif bad == 5:
            if prev == 2 and nxt != 4:
                result["G2"][depth] += 1
            elif prev == 1 and nxt != 4:
                result["G4"][depth] += 1
        elif bad == 0:
            prev2 = steps[t - 2] if t >= 2 else None
            if prev2 == 1 and prev == 2 and nxt != 1:
                result["D1"][depth] += 1
            elif prev2 == 1 and prev == 4 and nxt != 1:
                result["D2"][depth] += 1

    def visit(point, depth, a1, s1, a2, s2, a3, s3):
        result["counts"][depth] += 1
        if depth >= 2 and a1 != a2:
            result["turn"][depth] += 1
        if depth >= 3 and a1 == a3 and s1 == -s3 and a1 != a2:
            result["fold"][depth] += 1
        if classify and depth >= 1:
            classify_string(depth)
        if depth == max_depth:
            return
        for idx, delta in enumerate(PACKED_STEPS):
            nxt = point + delta
            if nxt in forbid or nxt in visited:
                continue
            visited.add(nxt)
            if classify:
                steps_hist.append(idx)
            visit(nxt, depth + 1, AXIS[idx], SIGN[idx], a1, s1, a2, s2)
            if classify:
                steps_hist.pop()
            visited.remove(nxt)

    visited = {start}
    visit(start, 0, -1, 0, -1, 0, -1, 0)
    return result


def count_height_restricted_walks(length: int, endpoint_height: int | None = None) -> int:
    """Independent tuple-coordinate enumerator for the finite slab families."""

    total = 0
    visited = {(0, 0, 0)}

    def visit(point: tuple[int, int, int], height: int, depth: int) -> None:
        nonlocal total
        if depth == length:
            total += 1
            return
        for dx, dy, dz in STEPS:
            nxt = (point[0] + dx, point[1] + dy, point[2] + dz)
            if nxt in visited:
                continue
            next_height = height - dx + dy + dz
            if endpoint_height is None:
                if next_height <= 0:
                    continue
            elif depth + 1 == length:
                if next_height != endpoint_height:
                    continue
            elif not 0 < next_height < endpoint_height:
                continue
            visited.add(nxt)
            visit(nxt, next_height, depth + 1)
            visited.remove(nxt)

    visit((0, 0, 0), 0, 0)
    return total


def walk_vertices(max_depth: int, forbidden: tuple[int, ...] = ()):
    out = [[] for _ in range(max_depth + 1)]
    forbid = frozenset(forbidden)
    start = pack(0, 0, 0)
    path = [start]

    def visit(point, depth):
        out[depth].append(tuple(path))
        if depth == max_depth:
            return
        for delta in PACKED_STEPS:
            nxt = point + delta
            if nxt in forbid or nxt in visited:
                continue
            visited.add(nxt)
            path.append(nxt)
            visit(nxt, depth + 1)
            path.pop()
            visited.remove(nxt)

    visited = {start}
    visit(start, 0)
    return out


def count_walks_from(start_pt, depth, forbidden_pts):
    if depth == 0:
        return 1
    total = 0
    visited = {start_pt}

    def visit(point, remaining):
        nonlocal total
        if remaining == 0:
            total += 1
            return
        for dx, dy, dz in STEPS:
            nxt = (point[0] + dx, point[1] + dy, point[2] + dz)
            if nxt in forbidden_pts or nxt in visited:
                continue
            visited.add(nxt)
            visit(nxt, remaining - 1)
            visited.remove(nxt)

    visit(start_pt, depth)
    return total


def defect_family_formulas(r: int) -> dict[str, int]:
    def suffix_count(j: int) -> int:
        return 1 if j == 0 else 2 * 3 ** (j - 1)

    g = sum(3 ** k * suffix_count(r - 2 - k) for k in range(0, r - 1)) if r >= 2 else 0
    d = sum(3 ** k * suffix_count(r - 3 - k) for k in range(0, r - 2)) if r >= 3 else 0
    return {"mono": 3 ** r, "G_each": g, "D_each": d, "L_r": 3 ** r + 4 * g + 2 * d}


def floor_decimal(value: str, places: int = 40) -> str:
    with localcontext() as context:
        context.prec = places + 20
        quantum = Decimal(1).scaleb(-places)
        return format(Decimal(value).quantize(quantum, rounding=ROUND_FLOOR), "f")


def atanh_neg_root_interval(magnitude: int, root: int):
    previous = mp.iv.dps
    try:
        mp.iv.dps = DPS
        x = mp.iv.mpf(magnitude)
        v = mp.iv.exp(mp.iv.log(x) * (mp.iv.mpf(-1) / root))
        text = str(mp.iv.log((1 + v) / (1 - v)) / 2).strip()
        lower, upper = text[1:-1].split(",", maxsplit=1)
        return lower.strip(), upper.strip()
    finally:
        mp.iv.dps = previous


def main() -> None:
    res_c = enumerate_walks(DEPTH)
    res_a = enumerate_walks(DEPTH, forbidden=(FORBID_A,), classify=True)
    res_b = enumerate_walks(DEPTH, forbidden=FORBID_B)
    c = res_c["counts"]
    t = res_c["turn"]
    d = res_c["fold"]
    a = res_a["counts"]
    aturn = res_a["turn"]
    b = res_b["counts"]

    check("c_series", tuple(c) == STORED_C, "independent backtracking reproduces c_0..c_10")
    check("a_series", tuple(a) == STORED_A, "independent backtracking reproduces a_0..a_10")
    check("b_series", tuple(b) == EXPECTED_B, "independent backtracking reproduces b_0..b_10")
    check("T_series", tuple(t) == EXPECTED_T, "independent turn counts T_0..T_10")
    check("D_series", tuple(d) == EXPECTED_D, "independent fold counts D_0..D_10")
    check("Aturn_series", tuple(aturn) == EXPECTED_ATURN, "independent avoid-and-turn counts")
    check(
        "turn_b_identity",
        all(t[m] == 24 * b[m - 2] for m in range(2, DEPTH + 1)),
        "time reversal gives T_m = 24 b_(m-2) for 2 <= m <= 10",
    )

    # These two counts are deliberately rebuilt from tuple coordinates, not
    # from the producer's packed walk enumerator.
    height_block_11_7 = count_height_restricted_walks(11, 7)
    halfspace_13 = count_height_restricted_walks(13)
    height_lower_35 = height_block_11_7 ** 2 * halfspace_13
    check(
        "height_slab_finite_counts",
        height_block_11_7 == HEIGHT_BLOCK_11_7 and halfspace_13 == HALFSPACE_13,
        "independent finite enumeration gives B_11,7=729000 and H_13=142016661",
    )
    check(
        "height_slab_a35_lower_bound",
        height_lower_35 == HEIGHT_LOWER_35 and height_lower_35 > defect_family_formulas(35)["L_r"],
        "two disjoint h-slabs plus a positive-h tail certify a_35 >= B^2 H",
    )

    source_hash = hashlib.sha256(SAW_SOURCE.read_bytes()).hexdigest()
    check("external_source_sha256", source_hash == C36_SHA256, "audited SBB2011 PDF hash")

    # ---- first-backtrack and multi-neighbour inequalities (m+n <= 10) ----
    ok_first = True
    ok_refined = True
    for total in range(2, DEPTH + 1):
        for n in range(1, total):
            m = total - n
            if not (c[total] <= c[m] * (c[n] - a[n - 1])):
                ok_first = False
            if m >= 2 and n >= 2:
                rhs = (
                    c[m] * (c[n] - a[n - 1])
                    - t[m] * b[n - 2]
                    - d[m] * (a[n - 1] - b[n - 2])
                )
                if not (c[total] <= rhs):
                    ok_refined = False
    check("first_backtrack_all_small", ok_first, "c_(m+n) <= c_m (c_n - a_(n-1)) for all m+n <= 10")
    check("multi_neighbour_all_small", ok_refined,
          "c_(m+n) <= c_m(c_n-a_(n-1)) - T_m b_(n-2) - D_m(a_(n-1)-b_(n-2)) for all m+n <= 10, m,n >= 2")

    # ---- exhaustive pair-level identity ----
    vert = walk_vertices(PAIR_TOTAL_MAX - 1)
    ok_pairs = True
    exact_splits = []
    for total in range(2, PAIR_TOTAL_MAX + 1):
        for n in range(1, total):
            m = total - n
            invalid = 0
            valid = 0
            for omega in vert[m]:
                wset = set(omega)
                shift = omega[-1]
                delta = pack_delta((shift >> 14) - OFF, ((shift >> 7) & 127) - OFF, (shift & 127) - OFF)
                for eta in vert[n]:
                    if any((p + delta) in wset for p in eta[1:]):
                        invalid += 1
                    else:
                        valid += 1
            fam = c[m] * a[n - 1]
            if n >= 2:
                fam += t[m] * b[n - 2] + d[m] * (a[n - 1] - b[n - 2])
            if valid != c[total] or invalid < fam:
                ok_pairs = False
            if invalid == fam and n >= 2:
                exact_splits.append((m, n))
    check("pair_level_bijection_and_cover", ok_pairs,
          "valid pairs biject with (m+n)-walks and invalid pairs cover all injected families for m+n <= 9")
    check("pair_level_exact_splits", set(exact_splits) >= {(2, 2), (3, 2)},
          f"families exhaust invalid pairs exactly at {sorted(set(exact_splits))}")

    # ---- triple inequalities (l+m+n <= 10) ----
    ok_triple = True
    beats_iteration = True
    for total in range(3, DEPTH + 1):
        for m in range(1, total - 1):
            for n in range(1, total - m):
                l = total - m - n
                if l < 1:
                    continue
                rhs7 = c[l] * (c[m] - a[m - 1]) * (c[n] - a[n - 1])
                if not (c[total] <= rhs7):
                    ok_triple = False
                if m >= 3 and n >= 2:
                    rhs8 = c[l] * (
                        (c[m] - a[m - 1]) * (c[n] - a[n - 1])
                        - (t[m] - aturn[m - 1]) * b[n - 2]
                    )
                    if not (c[total] <= rhs8) or not (rhs8 < rhs7):
                        ok_triple = False
                        beats_iteration = False
    check("triple_inequalities", ok_triple, "iterated and corner-corrected triple inequalities hold for l+m+n <= 10")
    check("triple_beats_iteration", beats_iteration,
          "the corner-corrected triple bound is strictly below pairwise iteration whenever applicable")

    # ---- chain instances (kn <= 10) ----
    ok_chain = True
    for n in range(1, 6):
        x_n = c[n] - a[n - 1]
        for k in range(2, DEPTH // n + 1):
            if c[k * n] > c[n] * x_n ** (k - 1):
                ok_chain = False
    check("chain_instances", ok_chain, "c_(kn) <= c_n (c_n - a_(n-1))^(k-1) for all kn <= 10")

    # ---- per-walk constancy (m <= 3, exhaustive) ----
    ok_perwalk = True
    tests = 0
    for m in range(2, 4):
        for omega in vert[m]:
            pts = [((p >> 14) - OFF, ((p >> 7) & 127) - OFF, (p & 127) - OFF) for p in omega]
            s1 = tuple(pts[m][i] - pts[m - 1][i] for i in range(3))
            s2 = tuple(pts[m - 1][i] - pts[m - 2][i] for i in range(3))
            turn = sum(x * y for x, y in zip(s1, s2)) == 0
            for n in range(2, 5):
                tests += 1
                if count_walks_from(pts[m - 1], n - 1, {pts[m]}) != a[n - 1]:
                    ok_perwalk = False
                if turn:
                    tests += 1
                    corner = tuple(pts[m][i] - s2[i] for i in range(3))
                    if count_walks_from(pts[m - 2], n - 2, {pts[m], corner}) != b[n - 2]:
                        ok_perwalk = False
    check("per_walk_universality", ok_perwalk, f"{tests} per-walk tail counts equal a/b exactly")

    # ---- defect families by brute force (r <= 10) ----
    ok_defect = True
    for r in range(0, DEPTH + 1):
        formulas = defect_family_formulas(r)
        g_counts = [res_a[name][r] for name in ("G1", "G2", "G3", "G4")]
        d_counts = [res_a[name][r] for name in ("D1", "D2")]
        union = res_a["mono"][r] + sum(g_counts) + sum(d_counts)
        if (
            res_a["mono"][r] != formulas["mono"]
            or any(gc != formulas["G_each"] for gc in g_counts)
            or any(dc != formulas["D_each"] for dc in d_counts)
            or union > a[r]
        ):
            ok_defect = False
    check("defect_families_bruteforce", ok_defect,
          "monotone and defect families match the exact position sums and stay below a_r for r <= 10")

    # ---- certificates from external data ----
    l_29 = defect_family_formulas(29)["L_r"]
    l_35 = defect_family_formulas(35)["L_r"]
    check("L_29_value", l_29 == 741377533262625, "L_29 = 3^29 + 4 G_29 + 2 D_29")
    check("L_35_value", l_35 == 644233352324156721, "L_35 = 3^35 + 4 G_35 + 2 D_35")

    best_base = EXTERNAL_COUNTS[36] - height_lower_35
    check("best_base_value", best_base == 2941295382858363225560670, "c_36 - H_35")
    check(
        "mu36_strictly_improves",
        best_base ** 36 < C36 ** 36,
        "(c_36 - H_35)^36 < c_36^36 as an integer comparison",
    )
    smaller_do_not_beat = all(
        not (EXTERNAL_COUNTS[n] - defect_family_formulas(n - 1)["L_r"]) ** 36 < C36 ** n
        for n in EXTERNAL_COUNTS
        if n < 36
    )
    check("smaller_external_do_not_beat", smaller_do_not_beat,
          "n = 30..35 chain instances do not beat c_36^(1/36) (finite-size deficit)")

    new_lo, _ = atanh_neg_root_interval(best_base, 36)
    old_lo, _ = atanh_neg_root_interval(C36, 36)
    new_floor = floor_decimal(new_lo)
    old_floor = floor_decimal(old_lo)
    check("incumbent_floor_reproduced", old_floor == INCUMBENT_K_LOWER,
          "atanh(c_36^(-1/36)) floors to the incumbent endpoint")
    check("new_floor_higher", Decimal(new_floor) > Decimal(old_floor),
          f"improved lower endpoint {new_floor}")

    # ---- n=30 lemma threshold ----
    needed = EXTERNAL_COUNTS[30] - EXTERNAL_COUNTS[32] // EXTERNAL_COUNTS[30]
    check("n30_threshold", needed == 270569905525454674592 and l_29 < needed,
          f"a_29 >= {needed} would prove mu^2 <= c_32/c_30; defect families give {l_29}")

    # ---- cross-check the stored artifact ----
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    data = artifact["data"]
    check(
        "artifact_agreement",
        data["main_certificate"]["best_base"] == best_base
        and data["main_certificate"]["best_lower_bound_name"] == "H_35"
        and data["main_certificate"]["best_lower_bound"] == height_lower_35
        and data["defect_family_lower_bound"]["L_29"] == l_29
        and data["defect_family_lower_bound"]["L_35"] == l_35
        and data["height_slab_lower_bound"]["H_35"] == height_lower_35
        and data["main_certificate"]["endpoint_comparison"]["new_lower_40f"] == new_floor
        and artifact["checks"] and all(row["passed"] for row in artifact["checks"]),
        "results/bounds/saw_refine.json matches this independent recomputation and all producer checks passed",
    )

    if FAILURES:
        print("\nFAILED:", FAILURES)
        raise SystemExit(1)
    print("OK")


if __name__ == "__main__":
    main()
