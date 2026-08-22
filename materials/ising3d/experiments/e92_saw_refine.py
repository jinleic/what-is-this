#!/usr/bin/env python3
"""Wave-10 first-backtrack SAW refinement certificate (e92).

Sharpen the first-backtrack submultiplicativity theorem of
``proofs/kc_interval2.md`` sec. 3 toward (but without assuming) the missing
finite-ratio lemma ``mu^2 <= c_32/c_30``.

Exact content of this script:

1. [THEOREM] multi-neighbour refinement of the injection.  After the forced
   first backtrack the m-walk still occupies a second *known* vertex near the
   junction (the square corner), and for fold endings even a third.  For
   m >= 2 and n >= 2,

       c_(m+n) <= c_m (c_n - a_(n-1)) - T_m b_(n-2) - D_m (a_(n-1) - b_(n-2))

   where T_m counts turn-ending m-walks, D_m counts fold-ending m-walks
   (last three steps s, t, -s with t perpendicular to s), and b_k counts
   k-walks avoiding the two vertices e_1 and e_1+e_2 forced by the corner
   geometry.  At (m,n) = (2,2) and (3,2) the injected families cover the
   invalid pairs *exactly* (equality holds).
2. [THEOREM] chained consequence: mu^n <= c_n - a_(n-1) for every n >= 1
   (induction on the all-m inequality plus the Fekete subsequence limit).
3. [THEOREM] explicit monotone/defect walk families giving a_r >= L_r with
   L_r = 3^r + 4 G_r + 2 D_r; G_r and D_r are exact position sums.
4. [THEOREM] a height-slab concatenation family gives a substantially stronger
   explicit lower bound on a_35 from two finite block enumerations.
5. [THEOREM] the chained n=36 instance, using the [EXTERNAL] exact count,
   strictly improves the incumbent mu <= c_36^(1/36) and hence the K_c lower
   endpoint.  Every comparison is a pure integer comparison.
6. [FALSIFIED] within the recorded direct-pair screening set, no instance
   beats c_36^(1/36); the triple inequality (for m >= 3) beats pairwise
   iteration in the checked finite cases but not numerically; the n=30 ratio
   lemma would follow from a_29 >= c_30 - floor(c_32/c_30), and the available
   defect-family bound is recorded against that threshold.

No decision uses a floating-point critical-coupling benchmark.  mpmath.iv
directed rounding at 90 dps formats endpoints only.
"""
from __future__ import annotations

from datetime import datetime, timezone
from decimal import Decimal, ROUND_FLOOR, localcontext
from fractions import Fraction
import hashlib
import json
import mpmath as mp
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))



SCRIPT = "experiments/e92_saw_refine.py"
RESULT_PATH = ROOT / "results" / "bounds" / "saw_refine.json"
SAW_SOURCE = ROOT / "sources" / "fulltext" / "schram_barkema_bisseling2011.pdf"
KC_BOUNDS_PATH = ROOT / "results" / "bounds" / "kc_bounds.json"
DPS = 90
REPORT_PLACES = 40
C36 = 2941370856334701726560670
C36_SHA256 = "898d8333e8d70acd279e29f226ff108550b132463aa662f86df590b0ad72cb12"
INCUMBENT_K_LOWER = "0.2122119011661678393310862783954278184914"
EXTERNAL_SAW_COUNTS = {
    30: 270569905525454674614,
    31: 1274191064726416905966,
    32: 5997359460809616886494,
    33: 28233744272563685150118,
    34: 132853629626823234210582,
    35: 625248129452557974777990,
    36: C36,
}
# In-repo stored exact series c_0..c_14 (e05 enumeration, cross-checked
# against OEIS A001412 upstream).  Independent anchor for our backtracking.
STORED_C_SMALL = (
    1, 6, 30, 150, 726, 3534, 16926, 81390, 387966, 1853886,
    8809878, 41934150, 198842742, 943971626, 4468911678,
)
STORED_A_SMALL = (1, 5, 25, 121, 589, 2821, 13565, 64661, 308981, 1468313)

MAX_C_DEPTH = 10
MAX_A_DEPTH = 10
MAX_B_DEPTH = 10
PAIR_TOTAL_MAX = 9          # exhaustive pair-level verification up to m+n = 9
PERWALK_M = 4               # per-walk symmetry verification up to m = 4
L_VERIFY_DEPTH = 10         # brute-force defect-family verification depth

# Fixed finite height-slab families used only for the a_35 certificate.
# The height is h(x,y,z)=-x+y+z.  A full block has all interior heights
# strictly between 0 and its endpoint height; a half-space tail has h>0.
HEIGHT_BLOCK_LENGTH = 11
HEIGHT_BLOCK_RISE = 7
HALFSPACE_TAIL_LENGTH = 13
HEIGHT_BLOCK_11_7 = 729000
HALFSPACE_13 = 142016661

# Packed lattice points: 7 bits per signed coordinate, offset +32.
OFF = 32
STEPS = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))


def pack(x: int, y: int, z: int) -> int:
    return ((x + OFF) << 14) | ((y + OFF) << 7) | (z + OFF)


def pack_delta(dx: int, dy: int, dz: int) -> int:
    # arithmetic sum: a bitwise OR is wrong when a lower field is negative
    return (dx << 14) + (dy << 7) + dz


PACKED_STEPS = tuple(pack_delta(*d) for d in STEPS)
AXIS = (0, 0, 1, 1, 2, 2)
SIGN = (1, -1, 1, -1, 1, -1)
FORBID_A = pack(1, 0, 0)                       # e_1
FORBID_B = (pack(1, 0, 0), pack(1, 1, 0))      # e_1 and e_1 + e_2
MONO_STEPS = (1, 2, 4)                          # -e_1, +e_2, +e_3

HEIGHT_DELTA = (-1, 1, 1, -1, 1, -1)            # h=-x+y+z under STEPS
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


def enumerate_walks(max_depth: int, forbidden: tuple[int, ...] = (), classify: bool = False) -> dict[str, list[int]]:
    """Exact DFS enumeration of rooted cubic walks avoiding `forbidden`.

    Returns per-depth totals, turn-ending and fold-ending counts, and with
    ``classify`` the defect-family class counts used for L_r.
    """
    result: dict[str, list[int]] = {
        key: [0] * (max_depth + 1) for key in ("counts", "turn", "fold")
    }
    for name in ("mono", "G1", "G2", "G3", "G4", "D1", "D2"):
        result[name] = [0] * (max_depth + 1)
    result["mono"][0] = 1  # the empty walk is monotone
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
        for s in steps[t + 1:]:
            if s not in MONO_STEPS:
                return  # more than one defect step: outside the families
        bad = steps[t]
        prev = steps[t - 1] if t >= 1 else None
        nxt = steps[t + 1] if t + 1 < depth else None
        if bad == 3:  # -e_2
            if prev == 4 and nxt != 2:
                result["G1"][depth] += 1
            elif prev == 1 and nxt != 2:
                result["G3"][depth] += 1
        elif bad == 5:  # -e_3
            if prev == 2 and nxt != 4:
                result["G2"][depth] += 1
            elif prev == 1 and nxt != 4:
                result["G4"][depth] += 1
        elif bad == 0:  # +e_1
            prev2 = steps[t - 2] if t >= 2 else None
            if prev2 == 1 and prev == 2 and nxt != 1:
                result["D1"][depth] += 1
            elif prev2 == 1 and prev == 4 and nxt != 1:
                result["D2"][depth] += 1

    def visit(point: int, depth: int, a1: int, s1: int, a2: int, s2: int, a3: int, s3: int) -> None:
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


def walk_vertices(max_depth: int, forbidden: tuple[int, ...] = ()) -> list[list[tuple[int, ...]]]:
    """All walks up to max_depth as tuples of packed vertices (start at 0)."""
    out: list[list[tuple[int, ...]]] = [[] for _ in range(max_depth + 1)]
    forbid = frozenset(forbidden)
    start = pack(0, 0, 0)
    path: list[int] = [start]

    def visit(point: int, depth: int) -> None:
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


def count_walks_from(start_pt: tuple[int, int, int], depth: int, forbidden_pts: set[tuple[int, int, int]]) -> int:
    """Exact number of `depth`-step walks from start_pt avoiding forbidden_pts."""
    if depth == 0:
        return 1
    total = 0
    visited = {start_pt}

    def visit(point: tuple[int, int, int], remaining: int) -> None:
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


def integer_strict_root(value: int, power: int) -> int:
    """Largest R with R**power < value (value > 1, power >= 1)."""
    if value <= 1:
        raise ValueError("value must exceed 1")
    low, high = 1, 1
    while high ** power < value:
        high *= 2
    while high - low > 1:
        mid = (low + high) // 2
        if mid ** power < value:
            low = mid
        else:
            high = mid
    return low


def defect_family_formulas(r: int) -> dict[str, int]:
    """Closed-form family counts L_r = 3^r + 4 G_r + 2 D_r."""

    def suffix_count(j: int) -> int:
        return 1 if j == 0 else 2 * 3 ** (j - 1)

    g = sum(3 ** k * suffix_count(r - 2 - k) for k in range(0, r - 1)) if r >= 2 else 0
    d = sum(3 ** k * suffix_count(r - 3 - k) for k in range(0, r - 2)) if r >= 3 else 0
    return {"mono": 3 ** r, "G_each": g, "D_each": d, "L_r": 3 ** r + 4 * g + 2 * d}


def count_height_restricted_walks(length: int, endpoint_height: int | None = None) -> int:
    """Count a finite self-avoiding height-slab or positive-half-space family.

    If ``endpoint_height`` is ``q``, every proper non-origin vertex has
    0 < h < q and the final vertex has h=q.  If it is ``None``, every
    non-origin vertex has h>0.  The coordinate packing remains safe at the
    fixed depths used here.
    """

    total = 0
    start = pack(0, 0, 0)
    visited = {start}

    def visit(point: int, height: int, depth: int) -> None:
        nonlocal total
        if depth == length:
            total += 1
            return
        for idx, delta in enumerate(PACKED_STEPS):
            nxt = point + delta
            if nxt in visited:
                continue
            next_height = height + HEIGHT_DELTA[idx]
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

    visit(start, 0, 0)
    return total


def main() -> None:
    started = time.monotonic()

    # ------------------------------------------------------------------
    # Enumeration of all exact series
    # ------------------------------------------------------------------
    enum_started = time.monotonic()
    res_c = enumerate_walks(MAX_C_DEPTH)
    res_a = enumerate_walks(MAX_A_DEPTH, forbidden=(FORBID_A,), classify=True)
    res_b = enumerate_walks(MAX_B_DEPTH, forbidden=FORBID_B)
    enumeration_seconds = time.monotonic() - enum_started

    height_enum_started = time.monotonic()
    height_block_11_7 = count_height_restricted_walks(HEIGHT_BLOCK_LENGTH, HEIGHT_BLOCK_RISE)
    halfspace_13 = count_height_restricted_walks(HALFSPACE_TAIL_LENGTH)
    height_enumeration_seconds = time.monotonic() - height_enum_started
    height_lower_35 = height_block_11_7 ** 2 * halfspace_13

    c_series = res_c["counts"]
    t_series = res_c["turn"]
    d_series = res_c["fold"]
    a_series = res_a["counts"]
    aturn_series = res_a["turn"]
    b_series = res_b["counts"]

    _check(
        "small_c_matches_stored_repo_series",
        tuple(c_series) == STORED_C_SMALL[: MAX_C_DEPTH + 1],
        "backtracking c_0..c_10 equals the OEIS-verified stored series",
    )
    _check(
        "small_a_matches_kc_interval2",
        tuple(a_series[: len(STORED_A_SMALL)]) == STORED_A_SMALL,
        "backtracking a_0..a_8 equals the kc_interval2 stored avoid-one-neighbour series",
    )
    stored_json = json.loads(KC_BOUNDS_PATH.read_text(encoding="utf-8"))
    stored_counts = stored_json["data"]["lower_bound_saw"]["enumeration"]["counts"]
    _check(
        "stored_c10_c14_json_agreement",
        all(int(stored_counts[str(i)]) == c_series[i] for i in range(MAX_C_DEPTH + 1)),
        "results/bounds/kc_bounds.json stores the same exact c_0..c_10",
    )
    source_hash = hashlib.sha256(SAW_SOURCE.read_bytes()).hexdigest()
    _check(
        "external_saw_source_hash",
        source_hash == C36_SHA256,
        "cached Schram--Barkema--Bisseling PDF has the audited SHA-256",
    )
    _check(
        "b_series_basic_sanity",
        b_series[0] == 1
        and b_series[1] == 5
        and b_series[2] == 24
        and all(b_series[r] <= a_series[r] for r in range(MAX_B_DEPTH + 1)),
        "b_0=1, b_1=5, b_2=24 by hand and b_r <= a_r throughout",
    )
    _check(
        "turn_b_identity",
        all(
            t_series[m] == 24 * b_series[m - 2]
            for m in range(2, min(MAX_C_DEPTH, MAX_B_DEPTH + 2) + 1)
        ),
        "time reversal and the 24 oriented turn prefixes give T_m = 24 b_(m-2) for 2 <= m <= 10",
    )
    _check(
        "height_slab_block_11_7",
        height_block_11_7 == HEIGHT_BLOCK_11_7,
        "direct finite enumeration of 11-step h-slab blocks ending at h=7",
    )
    _check(
        "height_halfspace_tail_13",
        halfspace_13 == HALFSPACE_13,
        "direct finite enumeration of 13-step strictly positive-h half-space walks",
    )
    _check(
        "height_slab_beats_defect_L35",
        height_lower_35 > defect_family_formulas(35)["L_r"],
        "two h-slab blocks followed by a positive-h tail give a_35 > L_35",
    )
    straight_series = [c_series[r] - t_series[r] if r >= 2 else 0 for r in range(MAX_C_DEPTH + 1)]
    _check(
        "straight_enders_inject_into_shorter_walks",
        all(straight_series[m] <= c_series[m - 1] for m in range(2, MAX_C_DEPTH + 1)),
        "S_m = c_m - T_m <= c_(m-1) by the prefix-injection for every m <= 10",
    )
    _check(
        "turn_beats_avoiding_turn",
        all(t_series[m] > aturn_series[m - 1] for m in range(3, 9)),
        "T_m > A^turn_(m-1) for 3 <= m <= 8, so the triple refinement strictly beats pairwise iteration",
    )

    counts: dict[int, int] = dict(enumerate(STORED_C_SMALL))
    counts.update(EXTERNAL_SAW_COUNTS)

    # ------------------------------------------------------------------
    # Theorem A: multi-neighbour refinement, exhaustive small verification
    # ------------------------------------------------------------------
    refined_rows: list[dict[str, object]] = []
    for total in range(4, MAX_C_DEPTH + 1):
        for n in range(2, total - 1):
            m = total - n
            if m < 2:
                continue
            lhs = counts[total]
            rhs = (
                counts[m] * (counts[n] - a_series[n - 1])
                - t_series[m] * b_series[n - 2]
                - d_series[m] * (a_series[n - 1] - b_series[n - 2])
            )
            passed = lhs <= rhs
            refined_rows.append(
                {"m": m, "n": n, "lhs_c_m_plus_n": lhs, "rhs": rhs, "passed": passed}
            )
            _check(
                f"refined_inequality_m{m}_n{n}",
                passed,
                "c_(m+n) <= c_m(c_n-a_(n-1)) - T_m b_(n-2) - D_m(a_(n-1)-b_(n-2))",
            )
    first_backtrack_rows: list[dict[str, object]] = []
    for total in range(2, MAX_C_DEPTH + 1):
        for n in range(1, total):
            m = total - n
            lhs = counts[total]
            rhs = counts[m] * (counts[n] - a_series[n - 1])
            first_backtrack_rows.append({"m": m, "n": n, "lhs": lhs, "rhs": rhs, "passed": lhs <= rhs})
            _check(f"first_backtrack_m{m}_n{n}", lhs <= rhs, "c_(m+n) <= c_m (c_n - a_(n-1))")

    # Triple inequalities: iterated pairwise and corner-corrected (m >= 3).
    triple_rows: list[dict[str, object]] = []
    for total in range(3, MAX_C_DEPTH + 1):
        for m in range(1, total - 1):
            for n in range(1, total - m):
                l = total - m - n
                if l < 1:
                    continue
                lhs = counts[total]
                rhs7 = counts[l] * (counts[m] - a_series[m - 1]) * (counts[n] - a_series[n - 1])
                _check(f"triple_iterated_l{l}_m{m}_n{n}", lhs <= rhs7,
                       "c_(l+m+n) <= c_l (c_m-a_(m-1))(c_n-a_(n-1))")
                if m >= 3 and n >= 2:
                    rhs8 = counts[l] * (
                        (counts[m] - a_series[m - 1]) * (counts[n] - a_series[n - 1])
                        - (t_series[m] - aturn_series[m - 1]) * b_series[n - 2]
                    )
                    _check(f"triple_refined_l{l}_m{m}_n{n}", lhs <= rhs8,
                           "triple with corner correction at the second junction")
                    triple_rows.append(
                        {"l": l, "m": m, "n": n, "lhs": lhs, "rhs_iterated": rhs7,
                         "rhs_triple": rhs8, "triple_strictly_beats_iteration": rhs8 < rhs7}
                    )
                else:
                    triple_rows.append({"l": l, "m": m, "n": n, "lhs": lhs, "rhs_iterated": rhs7})

    # ------------------------------------------------------------------
    # Exhaustive pair-level verification (m + n <= 9)
    # ------------------------------------------------------------------
    vert = walk_vertices(PAIR_TOTAL_MAX - 1)
    pair_rows: list[dict[str, object]] = []
    for total in range(2, PAIR_TOTAL_MAX + 1):
        for n in range(1, total):
            m = total - n
            invalid = 0
            valid = 0
            for omega in vert[m]:
                wset = set(omega)
                shift = omega[-1]
                dx = (shift >> 14) - OFF
                dy = ((shift >> 7) & 127) - OFF
                dz = (shift & 127) - OFF
                delta = pack_delta(dx, dy, dz)
                for eta in vert[n]:
                    if any((p + delta) in wset for p in eta[1:]):
                        invalid += 1
                    else:
                        valid += 1
            lhs = counts[total]
            fam_f1 = counts[m] * a_series[n - 1]
            fam_f2 = t_series[m] * b_series[n - 2] if n >= 2 else 0
            fam_f3 = d_series[m] * (a_series[n - 1] - b_series[n - 2]) if n >= 2 else 0
            injected = fam_f1 + fam_f2 + fam_f3
            pair_rows.append(
                {"m": m, "n": n, "c_m_c_n": counts[m] * counts[n], "c_m_plus_n": lhs,
                 "valid_pairs": valid, "invalid_pairs": invalid,
                 "family_f1": fam_f1, "family_f2": fam_f2, "family_f3": fam_f3,
                 "injected_total": injected, "slack": invalid - injected}
            )
            _check(f"pair_bijection_m{m}_n{n}", valid == lhs,
                   "valid concatenations biject with (m+n)-walks")
            _check(f"pair_injection_cover_m{m}_n{n}", invalid >= injected,
                   "invalid pairs cover all three injected families")
            if (m, n) in ((2, 2), (3, 2)):
                _check(
                    f"pair_injection_exact_m{m}_n{n}",
                    invalid == injected,
                    "at this split the three families exhaust the invalid pairs (equality)",
                )

    # ------------------------------------------------------------------
    # Per-walk symmetry verification (m <= 4)
    # ------------------------------------------------------------------
    perwalk_failures = 0
    perwalk_tests = 0
    for m in range(2, PERWALK_M + 1):
        for omega in vert[m]:
            pts = [
                ((p >> 14) - OFF, ((p >> 7) & 127) - OFF, (p & 127) - OFF)
                for p in omega
            ]
            s1 = tuple(pts[m][i] - pts[m - 1][i] for i in range(3))
            s2 = tuple(pts[m - 1][i] - pts[m - 2][i] for i in range(3))
            turn = sum(x * y for x, y in zip(s1, s2)) == 0
            for n in range(2, 6):
                got = count_walks_from(pts[m - 1], n - 1, {pts[m]})
                perwalk_tests += 1
                if got != a_series[n - 1]:
                    perwalk_failures += 1
                if turn:
                    corner = tuple(pts[m][i] - s2[i] for i in range(3))
                    got_b = count_walks_from(pts[m - 2], n - 2, {pts[m], corner})
                    perwalk_tests += 1
                    if got_b != b_series[n - 2]:
                        perwalk_failures += 1
                if m >= 3:
                    s3 = tuple(pts[m - 2][i] - pts[m - 3][i] for i in range(3))
                    if s3 == tuple(-x for x in s1) and turn:
                        got_f3 = count_walks_from(pts[m - 3], n - 1, {pts[m]})
                        perwalk_tests += 1
                        if got_f3 != a_series[n - 1]:
                            perwalk_failures += 1
                        got_ov = count_walks_from(pts[m - 2], n - 2, {pts[m - 3], pts[m]})
                        perwalk_tests += 1
                        if got_ov != b_series[n - 2]:
                            perwalk_failures += 1
    _check(
        "per_walk_tail_counts_universal",
        perwalk_failures == 0,
        f"all {perwalk_tests} per-walk tail counts equal a_(n-1)/b_(n-2) exactly (translation/rotation invariance)",
    )

    # ------------------------------------------------------------------
    # Defect-family lower bound a_r >= L_r: brute-force verification
    # ------------------------------------------------------------------
    defect_rows: list[dict[str, object]] = []
    for r in range(0, L_VERIFY_DEPTH + 1):
        formulas = defect_family_formulas(r)
        g_counts = [res_a[name][r] for name in ("G1", "G2", "G3", "G4")]
        d_counts = [res_a[name][r] for name in ("D1", "D2")]
        union = res_a["mono"][r] + sum(g_counts) + sum(d_counts)
        defect_rows.append(
            {"r": r, "a_r": a_series[r],
             "mono_bruteforce": res_a["mono"][r], "mono_formula": formulas["mono"],
             "G_each_bruteforce": g_counts, "G_each_formula": formulas["G_each"],
             "D_each_bruteforce": d_counts, "D_each_formula": formulas["D_each"],
             "counted_union_bruteforce": union, "L_r": formulas["L_r"]}
        )
        _check(f"defect_mono_r{r}", res_a["mono"][r] == formulas["mono"],
               "brute-force monotone-family count equals 3^r")
        _check(f"defect_G_r{r}", all(gc == formulas["G_each"] for gc in g_counts),
               "all four corner-defect families match the position sum")
        _check(f"defect_D_r{r}", all(dc == formulas["D_each"] for dc in d_counts),
               "both fold-defect families match the position sum")
        _check(f"defect_L_below_a_r{r}", union <= a_series[r],
               "counted union of defect families stays below a_r")

    # ------------------------------------------------------------------
    # Chained consequence and small chain instances
    # ------------------------------------------------------------------
    chain_rows: list[dict[str, object]] = []
    for n in range(1, 6):
        x_n = counts[n] - a_series[n - 1]
        for k in range(2, MAX_C_DEPTH // n + 1):
            if k * n > MAX_C_DEPTH:
                continue
            lhs = counts[k * n]
            rhs = counts[n] * x_n ** (k - 1)
            chain_rows.append({"n": n, "k": k, "lhs_c_kn": lhs, "rhs": rhs, "passed": lhs <= rhs})
            _check(f"chain_n{n}_k{k}", lhs <= rhs, "c_(kn) <= c_n (c_n - a_(n-1))^(k-1)")

    # ------------------------------------------------------------------
    # Main certificate: chain the first-backtrack theorem at external n.
    # The finite height-slab family is used at n=36; the defect family
    # supplies the comparison rows at the smaller stored external lengths.
    # ------------------------------------------------------------------
    external_lower_bounds = {
        n: defect_family_formulas(n - 1)["L_r"] for n in EXTERNAL_SAW_COUNTS
    }
    external_lower_names = {n: f"L_{n - 1}" for n in EXTERNAL_SAW_COUNTS}
    external_lower_bounds[36] = height_lower_35
    external_lower_names[36] = "H_35"
    external_rows: list[dict[str, object]] = []
    best_n: int | None = None
    best_base: int | None = None
    for n in sorted(EXTERNAL_SAW_COUNTS):
        lower = external_lower_bounds[n]
        base = EXTERNAL_SAW_COUNTS[n] - lower
        if base <= 0:
            raise AssertionError(f"nonpositive certificate base at n={n}")
        beats = base ** 36 < C36 ** n
        external_rows.append(
            {
                "n": n,
                "c_n": EXTERNAL_SAW_COUNTS[n],
                "lower_bound_name": external_lower_names[n],
                "lower_bound_on_a_(n-1)": lower,
                "mu_upper_base": base,
                "beats_c36_root_exactly": beats,
                "fingerprint_base36_minus_c36n": _fingerprint(base ** 36 - C36 ** n),
            }
        )
        if n == 36:
            _check(
                "external_chain_n36",
                beats,
                "mu^36 <= c_36 - H_35 strictly beats mu <= c_36^(1/36) as an integer comparison",
            )
        else:
            _check(
                f"external_chain_n{n}_does_not_beat",
                not beats,
                "at this smaller external n the finite-size deficit exceeds the defect-family correction (recorded)",
            )
        if best_base is None or base ** best_n < best_base ** n:
            best_n, best_base = n, base
    assert best_n is not None and best_base is not None
    _check(
        "best_external_instance_selection",
        all(
            not (EXTERNAL_SAW_COUNTS[n] - external_lower_bounds[n]) ** best_n
            < best_base ** n
            for n in EXTERNAL_SAW_COUNTS
        ),
        f"n = {best_n} attains the exact minimum of (c_n - recorded lower bound)^(1/n) over the stored external n",
    )

    new_lo, new_hi = _atanh_neg_root_interval(best_base, best_n)
    old_lo, old_hi = _atanh_neg_root_interval(C36, 36)
    new_floor = _floor_decimal(new_lo)
    old_floor = _floor_decimal(old_lo)
    _check(
        "incumbent_endpoint_reproduced",
        old_floor == INCUMBENT_K_LOWER,
        "atanh(c_36^(-1/36)) reproduces the incumbent lower endpoint at 40 places",
    )
    with localcontext() as context:
        context.prec = REPORT_PLACES + 20
        delta = Decimal(new_floor) - Decimal(old_floor)
    _check(
        "strict_endpoint_improvement",
        Decimal(new_floor) > Decimal(old_floor),
        f"new lower endpoint exceeds the incumbent by {delta}",
    )

    improvement_int = best_base ** 36 - C36 ** best_n

    # ------------------------------------------------------------------
    # Falsification records
    # ------------------------------------------------------------------
    direct_rows: list[dict[str, object]] = []
    best_direct: tuple[int, int, int, int] | None = None
    for m in list(range(2, 15)) + list(range(30, 37)):
        for n in range(1, 9):
            if m not in counts:
                continue
            correction = 0
            if n >= 2 and m <= MAX_C_DEPTH:
                correction = (
                    t_series[m] * b_series[n - 2]
                    + d_series[m] * (a_series[n - 1] - b_series[n - 2])
                )
            base_val = counts[m] * (counts[n] - a_series[n - 1]) - correction
            exponent = m + n
            beats = base_val ** 36 < C36 ** exponent
            needed = 0 if beats else base_val - integer_strict_root(C36 ** exponent, 36)
            direct_rows.append(
                {"m": m, "n": n, "base": base_val, "exponent": exponent,
                 "correction_applied": correction, "beats_c36_root": beats,
                 "shortfall_if_not": needed}
            )
            if not beats and (best_direct is None or needed < best_direct[0]):
                best_direct = (needed, m, n, base_val)
    assert best_direct is not None
    all_not_beat = all(not row["beats_c36_root"] for row in direct_rows)
    _check(
        "direct_pairs_do_not_beat_c36",
        all_not_beat,
        f"no direct pair instance beats c_36^(1/36); smallest shortfall {best_direct[0]} at (m,n)=({best_direct[1]},{best_direct[2]})",
    )

    l_ext, m_s, n_s = 30, 4, 2
    triple_base = counts[l_ext] * (
        (counts[m_s] - a_series[m_s - 1]) * (counts[n_s] - a_series[n_s - 1])
        - (t_series[m_s] - aturn_series[m_s - 1]) * b_series[n_s - 2]
    )
    chain_base = EXTERNAL_SAW_COUNTS[36] - external_lower_bounds[36]
    _check(
        "chain_beats_triple_at_certificate_scale",
        chain_base < triple_base,
        "the n=36 chained certificate is smaller than the best external-l triple certificate on mu^36",
    )

    c30, c32 = EXTERNAL_SAW_COUNTS[30], EXTERNAL_SAW_COUNTS[32]
    floor_ratio = c32 // c30
    needed_L = c30 - floor_ratio
    available_L = defect_family_formulas(29)["L_r"]
    lemma_gap = Fraction(needed_L, available_L)
    _check(
        "n30_lemma_gap_recorded",
        needed_L > available_L,
        f"proving a_29 >= c_30 - floor(c_32/c_30) = {needed_L} would give mu^2 <= c_32/c_30; defect families supply {available_L} (factor {float(lemma_gap):.1f} short)",
    )

    elapsed = time.monotonic() - started

    payload = {
        "provenance": {
            "script": SCRIPT,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": ".venv/bin/python",
            "arithmetic": "All decisions are Python integer comparisons; mpmath.iv with directed rounding at 90 dps formats endpoints only.",
            "mpmath_iv_dps": DPS,
            "elapsed_seconds": elapsed,
            "enumeration_seconds": enumeration_seconds,
            "height_enumeration_seconds": height_enumeration_seconds,
        },
        "data": {
            "classification": (
                "[THEOREM] the first-backtrack chain with the explicit height-slab lower bound "
                "on a_35 strictly improves the connective-constant upper bound and the K_c lower "
                "endpoint; the multi-neighbour corner/fold corrections are exact but a direct pair "
                "instance still cannot beat c_36^(1/36); mu^2 <= c_32/c_30 remains open with the "
                "exact deficit recorded."
            ),
            "series": {
                "c_0_to_10": c_series,
                "a_0_to_10": a_series,
                "b_0_to_10": b_series,
                "T_0_to_10": t_series,
                "D_0_to_10": d_series,
                "Aturn_0_to_10": aturn_series,
                "external_c_30_to_36": {str(k): v for k, v in EXTERNAL_SAW_COUNTS.items()},
            },
            "source": {
                "path": "sources/fulltext/schram_barkema_bisseling2011.pdf",
                "sha256": source_hash,
                "citation": "Schram, Barkema, Bisseling (2011), J. Stat. Mech. P06019, Table I",
            },
            "theorem_multi_neighbour": {
                "status": "[THEOREM]",
                "statement": (
                    "For m >= 2 and n >= 2: c_(m+n) <= c_m (c_n - a_(n-1)) - T_m b_(n-2) "
                    "- D_m (a_(n-1) - b_(n-2)); T_m counts m-walks whose last two steps are "
                    "perpendicular, D_m those whose last three steps are s, t, -s with t "
                    "perpendicular to s, and b_k counts k-walks avoiding e_1 and e_1+e_2."
                ),
                "injections": [
                    "F1 (first backtrack, prior wave): eta_1 = omega_(m-1); tail from omega_(m-1) avoiding omega_m: a_(n-1) per omega.",
                    "F2 (square corner, NEW): eta_1 = omega_m - s_2, eta_2 = omega_(m-2); tail from omega_(m-2) avoiding {omega_m, omega_m - s_2}: b_(n-2) per turn-ending omega.",
                    "F3 (fold return, NEW): eta_1 = omega_(m-3) (fold endings only); tail from omega_(m-3) avoiding omega_m: a_(n-1) per fold-ending omega, of which b_(n-2) are already in F2, leaving a_(n-1) - b_(n-2).",
                ],
                "exact_splits": [
                    {"m": row["m"], "n": row["n"]}
                    for row in pair_rows
                    if row["slack"] == 0 and row["n"] >= 2
                ],
                "small_checks": refined_rows,
                "first_backtrack_checks": first_backtrack_rows,
                "pair_level_exhaustive": pair_rows,
                "triple_checks": triple_rows,
            },
            "lemma_turn_b_identity": {
                "status": "[LEMMA]",
                "statement": "For every m >= 2, T_m = 24 b_(m-2).",
                "proof_sketch": (
                    "Reverse a turn-ending walk. Its first two directed steps are one of 24 ordered "
                    "perpendicular pairs; after translation and a cubic-lattice symmetry the remaining "
                    "tail avoids the fixed corner pair e_1,e_1+e_2 and is counted by b_(m-2)."
                ),
                "finite_checks": {"range": "2 <= m <= 10", "values": t_series[2:]},
            },
            "theorem_chain": {
                "status": "[THEOREM]",
                "statement": "For every n >= 1: mu^n <= c_n - a_(n-1).",
                "proof_sketch": (
                    "Theorem (4) of kc_interval2 holds for every m at fixed n; induction gives "
                    "c_(kn) <= c_n (c_n - a_(n-1))^(k-1); the Fekete subsequence limit along "
                    "multiples of n equals mu."
                ),
                "small_chain_checks": chain_rows,
            },
            "defect_family_lower_bound": {
                "status": "[THEOREM]",
                "statement": (
                    "a_r >= L_r := 3^r + 4 G_r + 2 D_r where G_r = sum_k 3^k s(r-2-k), "
                    "D_r = sum_k 3^k s(r-3-k), s(0)=1, s(j)=2*3^(j-1); the six families are "
                    "pairwise disjoint monotone-step walks with at most one defect step, the "
                    "defect being -e_2/-e_3 preceded by e_3/e_2/-e_1 or +e_1 preceded by "
                    "(-e_1, e_2/e_3), with the suffix's first step restricted away from one "
                    "direction."
                ),
                "bruteforce_rows": defect_rows,
                "L_29": defect_family_formulas(29)["L_r"],
                "L_35": defect_family_formulas(35)["L_r"],
            },
            "height_slab_lower_bound": {
                "status": "[THEOREM]",
                "finite_count_status": "[COMPUTATION]",
                "statement": (
                    "With h(x,y,z)=-x+y+z, let B be the 11-step self-avoiding blocks from h=0 "
                    "to h=7 with all interior heights strictly between, and let H be the 13-step "
                    "self-avoiding walks whose non-origin heights are positive. Then a_35 >= B^2 H."
                ),
                "block_count_B_11_7": height_block_11_7,
                "tail_count_H_13": halfspace_13,
                "H_35": height_lower_35,
                "defect_L_35_for_comparison": defect_family_formulas(35)["L_r"],
                "construction": (
                    "Concatenate two translated B blocks at heights 0..7 and 7..14, then a "
                    "translated H tail above height 14. Their interiors are disjoint and meet only "
                    "at prescribed splice endpoints; h(e_1)=-1 proves avoidance of e_1."
                ),
            },
            "main_certificate": {
                "status": "[THEOREM]",
                "external_input_status": "[EXTERNAL]",
                "statement": (
                    "mu^n <= c_n - L_(n-1) for stored n=30,...,35 and "
                    "mu^36 <= c_36 - H_35; the recorded n=36 instance improves "
                    "mu <= c_36^(1/36)."
                ),
                "rows": external_rows,
                "best_n": best_n,
                "best_base": best_base,
                "best_lower_bound_name": external_lower_names[best_n],
                "best_lower_bound": external_lower_bounds[best_n],
                "improvement_integer": {
                    "expression": f"(c_{best_n} - {external_lower_names[best_n]})^(36) - (c_36)^({best_n})",
                    "fingerprint": _fingerprint(improvement_int),
                },
                "endpoint_comparison": {
                    "old_lower_40f": old_floor,
                    "new_lower_40f": new_floor,
                    "new_interval_90dps": [new_lo, new_hi],
                    "old_interval_90dps": [old_lo, old_hi],
                    "strictly_improves": True,
                },
            },
            "falsifications": {
                "direct_pairs": {
                    "status": "[FALSIFIED]",
                    "scope": "as an improvement route at enumerable scale",
                    "rows": direct_rows,
                    "smallest_shortfall": {
                        "needed_reduction": best_direct[0],
                        "at_m": best_direct[1],
                        "at_n": best_direct[2],
                        "base": best_direct[3],
                    },
                },
                "triple_vs_chain": {
                    "status": "[COMPUTATION]",
                    "triple_instance": {"l": l_ext, "m": m_s, "n": n_s, "mu^36_upper": triple_base},
                    "chain_instance": {"n": 36, "mu^36_upper": chain_base},
                    "finding": "the triple inequality strictly improves pairwise iteration in the enumerated m>=3 cases but is far weaker than the chained certificate",
                },
                "missing_n30_ratio_lemma": {
                    "status": "[UNRESOLVED]",
                    "sufficient_threshold": needed_L,
                    "available_defect_bound": available_L,
                    "gap_factor_fraction": f"{lemma_gap.numerator}/{lemma_gap.denominator}",
                    "statement": "a_29 >= c_30 - floor(c_32/c_30) together with mu^30 <= c_30 - a_29 would prove mu^2 <= c_32/c_30.",
                },
            },
        },
        "checks": checks,
    }

    RESULT_PATH.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    print(f"\nwrote {RESULT_PATH}")
    print(f"elapsed {elapsed:.1f}s")
    print(f"new lower endpoint (floor 40): {new_floor}")
    print(f"old lower endpoint (floor 40): {old_floor}")


if __name__ == "__main__":
    try:
        main()
    except AssertionError as exc:
        print(f"CERTIFICATE FAILURE: {exc}")
        raise
