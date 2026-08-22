#!/usr/bin/env python3
"""Clean-room standalone verifier for ``experiments/e116_saw_union``.

Independently recomputes every constant of the wave-11 union certificate
with structurally different implementations where practical:

* iterative explicit-stack DFS on tuple coordinates (the producer uses
  recursion on packed integers);
* the coordinate-monotone count by a closed-form per-pattern
  inclusion-exclusion over surjective words (the producer uses a
  27-state automaton);
* the C cap H overlap by per-pattern Moebius inversion whose per-orthant
  block/tail counts come from restricted visited-set enumerations (the
  producer aggregates globally and counts by height DP);
* the endpoint by a 150-dps enclosure (the producer uses 120 dps).

It also verifies set overlaps/disjointness, the integer identities, the
submultiplicative chain metadata, and the directed 40-place floor, and it
fails on plausible off-by-one/overlap defects (wrong block interiors,
forgotten overlap subtraction, double counting, mis-rounded floors).

Run:  PYTHONPATH=src .venv/bin/python tests/test_saw_union.py
"""
from __future__ import annotations

from decimal import Decimal, ROUND_FLOOR, localcontext
from itertools import product
import json
from pathlib import Path

import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "bounds" / "saw_union.json"
SAW_SOURCE = ROOT / "sources" / "fulltext" / "schram_barkema_bisseling2011.pdf"

C36 = 2941370856334701726560670
C36_SHA256 = "898d8333e8d70acd279e29f226ff108550b132463aa662f86df590b0ad72cb12"
DPS = 150
SERIES_DEPTH = 9
CM_BRUTE_DEPTH = 7
L_BRUTE_DEPTH = 9
BLOCK_LEN, BLOCK_RISE, TAIL_LEN = 11, 7, 13
MINI_CONFIGS = ((5, 3, 4), (5, 3, 5), (5, 3, 6))

STEPS = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))
DH = (-1, 1, 1, -1, 1, -1)
MONO = (1, 2, 4)

FAILURES: list[str] = []


def check(name: str, passed: bool, detail: str = "") -> None:
    print(("PASS" if passed else "FAIL") + f": {name}" + (f" ({detail})" if detail else ""))
    if not passed:
        FAILURES.append(name)


# ----------------------------------------------------------------------
# Iterative tuple-coordinate enumerators.
# ----------------------------------------------------------------------



def iter_walk_words(max_depth: int) -> list[list[tuple[int, ...]]]:
    """e_1-avoiding SAW step-index words per length, iterative DFS."""
    out: list[list[tuple[int, ...]]] = [[] for _ in range(max_depth + 1)]
    out[0].append(())
    start = (0, 0, 0)
    stack = [(start, (start,), ())]
    while stack:
        point, vertices, word = stack.pop()
        depth = len(word)
        if depth >= 1:
            out[depth].append(word)
        if depth == max_depth:
            continue
        x, y, z = point
        for i, (dx, dy, dz) in enumerate(STEPS):
            q = (x + dx, y + dy, z + dz)
            if q == (1, 0, 0) or q in vertices:
                continue
            stack.append((q, vertices + (q,), word + (i,)))
    return out


def iter_height_block(length: int, rise: int) -> int:
    """11-or-fewer-step blocks: `length` steps, h 0->rise, interior (0,rise)."""
    total = 0
    start = (0, 0, 0)
    stack = [(start, 0, 0, (start,))]
    while stack:
        point, height, depth, vertices = stack.pop()
        if depth == length:
            total += 1
            continue
        x, y, z = point
        for i, (dx, dy, dz) in enumerate(STEPS):
            q = (x + dx, y + dy, z + dz)
            if q in vertices:
                continue
            nh = height + DH[i]
            if depth + 1 == length:
                if nh != rise:
                    continue
            elif not 0 < nh < rise:
                continue
            stack.append((q, nh, depth + 1, vertices + (q,)))
    return total


def iter_halfspace_tail(length: int) -> int:
    total = 0
    start = (0, 0, 0)
    stack = [(start, 0, 0, (start,))]
    while stack:
        point, height, depth, vertices = stack.pop()
        if depth == length:
            total += 1
            continue
        x, y, z = point
        for i, (dx, dy, dz) in enumerate(STEPS):
            q = (x + dx, y + dy, z + dz)
            if q in vertices:
                continue
            if height + DH[i] <= 0:
                continue
            stack.append((q, height + DH[i], depth + 1, vertices + (q,)))
    return total


def iter_restricted_block(length: int, rise: int, letters: tuple[int, ...]) -> int:
    allowed = frozenset(letters)
    total = 0
    start = (0, 0, 0)
    stack = [(start, 0, 0, (start,))]
    while stack:
        point, height, depth, vertices = stack.pop()
        if depth == length:
            total += 1
            continue
        x, y, z = point
        for i in allowed:
            dx, dy, dz = STEPS[i]
            q = (x + dx, y + dy, z + dz)
            if q in vertices:
                continue
            nh = height + DH[i]
            if depth + 1 == length:
                if nh != rise:
                    continue
            elif not 0 < nh < rise:
                continue
            stack.append((q, nh, depth + 1, vertices + (q,)))
    return total


def iter_restricted_tail(length: int, letters: tuple[int, ...]) -> int:
    allowed = frozenset(letters)
    total = 0
    start = (0, 0, 0)
    stack = [(start, 0, 0, (start,))]
    while stack:
        point, height, depth, vertices = stack.pop()
        if depth == length:
            total += 1
            continue
        x, y, z = point
        for i in allowed:
            dx, dy, dz = STEPS[i]
            q = (x + dx, y + dy, z + dz)
            if q in vertices:
                continue
            if height + DH[i] <= 0:
                continue
            stack.append((q, height + DH[i], depth + 1, vertices + (q,)))
    return total


# ----------------------------------------------------------------------
# Family predicates.
# ----------------------------------------------------------------------

def is_cm_word(word: tuple[int, ...]) -> bool:
    slots = [0, 0, 0]
    for i in word:
        ax = i // 2
        slot = 2 if i % 2 == 0 else 1
        if slots[ax] != 0 and slots[ax] != slot:
            return False
        slots[ax] = slot
    return True


def classify_l(word: tuple[int, ...]) -> str | None:
    for s in word:
        if s not in MONO:
            break
    else:
        return "mono"
    t = next(i for i, s in enumerate(word) if s not in MONO)
    for s in word[t + 1:]:
        if s not in MONO:
            return None
    bad = word[t]
    prev = word[t - 1] if t >= 1 else None
    nxt = word[t + 1] if t + 1 < len(word) else None
    if bad == 3:
        if prev == 4 and nxt != 2:
            return "G1"
        if prev == 1 and nxt != 2:
            return "G3"
    elif bad == 5:
        if prev == 2 and nxt != 4:
            return "G2"
        if prev == 1 and nxt != 4:
            return "G4"
    elif bad == 0:
        prev2 = word[t - 2] if t >= 2 else None
        if prev2 == 1 and prev == 2 and nxt != 1:
            return "D1"
        if prev2 == 1 and prev == 4 and nxt != 1:
            return "D2"
    return None


def defect_L(r: int) -> int:
    def s(j: int) -> int:
        return 1 if j == 0 else 2 * 3 ** (j - 1)

    g = sum(3 ** k * s(r - 2 - k) for k in range(0, r - 1)) if r >= 2 else 0
    d = sum(3 ** k * s(r - 3 - k) for k in range(0, r - 2)) if r >= 3 else 0
    return 3 ** r + 4 * g + 2 * d


# ----------------------------------------------------------------------
# Closed-form coordinate-monotone count (independent of the automaton).
# ----------------------------------------------------------------------

def cm_closed_form(r: int) -> int:
    """|C(r)| by per-pattern surjective-word inclusion-exclusion."""
    if r == 0:
        return 1  # the unique empty word has the empty used-letter pattern
    total = 0
    for code in product((None, 0, 1), repeat=3):
        letters = []
        for ax, s in enumerate(code):
            if s is None:
                continue
            letters.append(ax * 2 + (0 if s == 0 else 1))
        if not letters:
            continue
        P = tuple(letters)
        # words over P using every letter of P
        all_used = 0
        for mask in range(1 << len(P)):
            sub = [P[j] for j in range(len(P)) if mask >> j & 1]
            all_used += (-1) ** (len(P) - len(sub)) * (len(sub) ** r if sub else 0)
        # ... and with first letter +e_1 when +e_1 in P.  The suffix must
        # contain every other letter, but +e_1 remains available there.
        first_plus = 0
        if 0 in P:
            rest = [i for i in P if i != 0]
            for mask in range(1 << len(rest)):
                sub = [rest[j] for j in range(len(rest)) if mask >> j & 1]
                first_plus += (-1) ** (len(rest) - len(sub)) * (1 + len(sub)) ** (r - 1)
        total += all_used - first_plus
    return total


# ----------------------------------------------------------------------
# Mini slab families, iterative.
# ----------------------------------------------------------------------

def mini_slab_cm_count(blk_len: int, blk_rise: int, tail_len: int) -> tuple[int, int]:
    """Returns the full mini-family count and its coordinate-monotone members."""
    r = 2 * blk_len + tail_len
    total_family = 0
    total_cm = 0
    start = (0, 0, 0)

    def allowed(nh: int, depth_next: int) -> bool:
        if depth_next == blk_len:
            return nh == blk_rise
        if depth_next < blk_len:
            return 0 < nh < blk_rise
        if depth_next == 2 * blk_len:
            return nh == 2 * blk_rise
        if depth_next < 2 * blk_len:
            return blk_rise < nh < 2 * blk_rise
        return nh > 2 * blk_rise

    stack = [(start, 0, 0, (start,), (0, 0, 0), True)]
    while stack:
        point, height, depth, vertices, slots, cm_ok = stack.pop()
        if depth == r:
            total_family += 1
            if cm_ok:
                total_cm += 1
            continue
        x, y, z = point
        for i, (dx, dy, dz) in enumerate(STEPS):
            q = (x + dx, y + dy, z + dz)
            if q in vertices:
                continue
            nh = height + DH[i]
            if not allowed(nh, depth + 1):
                continue
            ax = i // 2
            slot = 2 if i % 2 == 0 else 1
            ok = cm_ok and (slots[ax] == 0 or slots[ax] == slot)
            ns = list(slots)
            ns[ax] = slot
            stack.append((q, nh, depth + 1, vertices + (q,), tuple(ns), ok))
    return total_family, total_cm

def mini_slab_cm_pruned(blk_len: int, blk_rise: int, tail_len: int) -> int:
    """Direct coordinate-monotone mini-slab count, pruning sign conflicts."""
    r = 2 * blk_len + tail_len
    total = 0
    start = (0, 0, 0)

    def allowed(nh: int, depth_next: int) -> bool:
        if depth_next == blk_len:
            return nh == blk_rise
        if depth_next < blk_len:
            return 0 < nh < blk_rise
        if depth_next == 2 * blk_len:
            return nh == 2 * blk_rise
        if depth_next < 2 * blk_len:
            return blk_rise < nh < 2 * blk_rise
        return nh > 2 * blk_rise

    stack = [(start, 0, 0, (start,), (0, 0, 0))]
    while stack:
        point, height, depth, vertices, slots = stack.pop()
        if depth == r:
            total += 1
            continue
        x, y, z = point
        for i, (dx, dy, dz) in enumerate(STEPS):
            q = (x + dx, y + dy, z + dz)
            if q in vertices:
                continue
            nh = height + DH[i]
            if not allowed(nh, depth + 1):
                continue
            ax = i // 2
            slot = 2 if i % 2 == 0 else 1
            if slots[ax] != 0 and slots[ax] != slot:
                continue
            ns = list(slots)
            ns[ax] = slot
            stack.append((q, nh, depth + 1, vertices + (q,), tuple(ns)))
    return total


def mini_slab_h_cap_l(blk_len: int, blk_rise: int, tail_len: int) -> int:
    r = 2 * blk_len + tail_len
    total = 0
    start = (0, 0, 0)

    def allowed(nh: int, depth_next: int) -> bool:
        if depth_next == blk_len:
            return nh == blk_rise
        if depth_next < blk_len:
            return 0 < nh < blk_rise
        if depth_next == 2 * blk_len:
            return nh == 2 * blk_rise
        if depth_next < 2 * blk_len:
            return blk_rise < nh < 2 * blk_rise
        return nh > 2 * blk_rise

    stack = [(start, 0, 0, (start,), ())]
    while stack:
        point, height, depth, vertices, word = stack.pop()
        if depth == r:
            if classify_l(word) is not None:
                total += 1
            continue
        x, y, z = point
        for i, (dx, dy, dz) in enumerate(STEPS):
            q = (x + dx, y + dy, z + dz)
            if q in vertices:
                continue
            nh = height + DH[i]
            if not allowed(nh, depth + 1):
                continue
            if sum(1 for s in word + (i,) if s not in MONO) > 1:
                continue  # any L member has at most one defect letter
            stack.append((q, nh, depth + 1, vertices + (q,), word + (i,)))
    return total


# ----------------------------------------------------------------------
# Per-pattern Moebius for C cap H (independent aggregation).
# ----------------------------------------------------------------------

def c_cap_h_per_pattern(blk_len: int, blk_rise: int, tail_len: int) -> int:
    total = 0
    for code in product((None, 0, 1), repeat=3):
        letters = []
        for ax, s in enumerate(code):
            if s is None:
                continue
            letters.append(ax * 2 + (0 if s == 0 else 1))
        if not letters:
            continue
        P = tuple(letters)
        n_p = 0
        for mask in range(1 << len(P)):
            sub = tuple(P[j] for j in range(len(P)) if mask >> j & 1)
            if not sub:
                continue
            b1 = iter_restricted_block(blk_len, blk_rise, sub)
            tq = iter_restricted_tail(tail_len, sub)
            n_p += (-1) ** (len(P) - len(sub)) * b1 * b1 * tq
        total += n_p
    return total


# ----------------------------------------------------------------------
# Interval helpers.
# ----------------------------------------------------------------------

def enclosure(magnitude: int, root: int) -> tuple[str, str]:
    with mp.workdps(DPS):
        mp.iv.dps = DPS
        x = mp.iv.mpf(magnitude)
        v = mp.iv.exp(mp.iv.log(x) * (mp.iv.mpf(-1) / root))
        text = str(mp.iv.log((1 + v) / (1 - v)) / 2).strip()
        lo, hi = text[1:-1].split(",")
        return lo.strip(), hi.strip()


def floor40(value: str) -> str:
    with localcontext() as ctx:
        ctx.prec = 80
        return format(Decimal(value).quantize(Decimal(1).scaleb(-40), rounding=ROUND_FLOOR), "f")


def main() -> None:
    import hashlib

    # ------------------------------------------------------------------
    # External source audit
    # ------------------------------------------------------------------
    check("source_hash", hashlib.sha256(SAW_SOURCE.read_bytes()).hexdigest() == C36_SHA256,
          "cached SAW source PDF hash")

    # ------------------------------------------------------------------
    # Series by independent iterative enumeration
    # ------------------------------------------------------------------
    words = iter_walk_words(SERIES_DEPTH)
    # c_r needs ALL SAWs (not avoiding e_1): recompute without the e_1 filter
    c_all = [0] * (SERIES_DEPTH + 1)
    stack = [((0, 0, 0), 0, ((0, 0, 0),))]
    while stack:
        point, depth, vertices = stack.pop()
        c_all[depth] += 1
        if depth == SERIES_DEPTH:
            continue
        x, y, z = point
        for dx, dy, dz in STEPS:
            q = (x + dx, y + dy, z + dz)
            if q not in vertices:
                stack.append((q, depth + 1, vertices + (q,)))
    a_series = [len(words[r]) for r in range(SERIES_DEPTH + 1)]
    check("c_series", tuple(c_all) == (1, 6, 30, 150, 726, 3534, 16926, 81390, 387966, 1853886),
          "iterative c_0..c_9 enumeration")
    check("a_series", tuple(a_series[:10]) == (1, 5, 25, 121, 589, 2821, 13565, 64661, 308981, 1468313),
          "iterative a_0..a_9 enumeration")

    # ------------------------------------------------------------------
    # Family H
    # ------------------------------------------------------------------
    B = iter_height_block(BLOCK_LEN, BLOCK_RISE)
    T = iter_halfspace_tail(TAIL_LEN)
    H35 = B * B * T
    check("B_11_7_iterative", B == 729000, f"blocks: {B}")
    check("T_13_iterative", T == 142016661, f"tails: {T}")
    check("H35_iterative", H35 == 75473476338501000000, "H_35 = B^2 T")

    # ------------------------------------------------------------------
    # Family L and small-scale brute-force agreements
    # ------------------------------------------------------------------
    L35 = defect_L(35)
    check("L35_independent", L35 == 644233352324156721, "defect family L_35")
    ok_rows = True
    cm_cap_l_small: dict[int, int] = {}
    for r in range(L_BRUTE_DEPTH + 1):
        fam: dict[str, int] = {}
        cmf: dict[str, int] = {}
        for w in words[r]:
            f = classify_l(w)
            if f is not None:
                fam[f] = fam.get(f, 0) + 1
                if is_cm_word(w):
                    cmf[f] = cmf.get(f, 0) + 1
        def s_(j: int) -> int:
            return 1 if j == 0 else 2 * 3 ** (j - 1)
        g = sum(3 ** k * s_(r - 2 - k) for k in range(0, r - 1)) if r >= 2 else 0
        d = sum(3 ** k * s_(r - 3 - k) for k in range(0, r - 2)) if r >= 3 else 0
        exp_cm_g = (r - 1) * 2 ** (r - 2) if r >= 2 else 0
        if not (fam.get("mono", 0) == 3 ** r
                and all(fam.get(f"G{i}", 0) == g for i in (1, 2, 3, 4))
                and all(fam.get(f"D{i}", 0) == d for i in (1, 2))
                and cmf.get("mono", 0) == 3 ** r
                and all(cmf.get(f"G{i}", 0) == exp_cm_g for i in (1, 2, 3, 4))
                and all(cmf.get(f"D{i}", 0) == 0 for i in (1, 2))):
            ok_rows = False
        cm_cap_l_small[r] = sum(cmf.values())
    check("L_families_small_r", ok_rows, "family and CM-overlap counts r<=9 by classification")
    check("cm_cap_L_small_formula",
          all(cm_cap_l_small[r] == 3 ** r + 4 * ((r - 1) * 2 ** (r - 2) if r >= 2 else 0)
              for r in range(L_BRUTE_DEPTH + 1)),
          "C cap L (r) = 3^r + 4 (r-1) 2^(r-2) for r<=9")

    # ------------------------------------------------------------------
    # Family C
    # ------------------------------------------------------------------
    ok_cm = True
    for r in range(CM_BRUTE_DEPTH + 1):
        brute = sum(1 for w in words[r] if is_cm_word(w))
        if brute != cm_closed_form(r):
            ok_cm = False
    check("cm_closed_form_small_r", ok_cm,
          "closed-form CM count matches exhaustive classification r<=7")
    C35 = cm_closed_form(35)
    check("C35_closed_form", C35 == 333543290395947705, f"|C| = {C35}")

    # ------------------------------------------------------------------
    # Overlaps
    # ------------------------------------------------------------------
    c_cap_l = 3 ** 35 + 4 * sum(2 ** 33 for _ in range(1, 35))
    check("C_cap_L_independent", c_cap_l == 50032713330104219,
          "C cap L = 3^35 + 4 * sum_{t=1..34} 2^33")
    c_cap_h = c_cap_h_per_pattern(BLOCK_LEN, BLOCK_RISE, TAIL_LEN)
    check("C_cap_H_per_pattern", c_cap_h == 157640944278888,
          f"C cap H by per-pattern Moebius + restricted DFS = {c_cap_h}")

    for (bl, br, tl) in MINI_CONFIGS:
        cm_n = mini_slab_cm_pruned(bl, br, tl)
        mo = c_cap_h_per_pattern(bl, br, tl)
        hl = mini_slab_h_cap_l(bl, br, tl)
        detail = f"mini ({bl},{br},{tl}): direct CM={cm_n}, Moebius={mo}"
        if (bl, br, tl) == (5, 3, 4):
            fam_n, full_cm_n = mini_slab_cm_count(bl, br, tl)
            check("mini_534_full_vs_pruned", full_cm_n == cm_n,
                  "full mini-family traversal and sign-pruned traversal agree")
            check("mini_534_size", fam_n == 108 * 108 * 171,
                  "|H(5,3,4)| = B(5,3)^2 T(4)")
            detail = f"{detail}, |H|={fam_n}"
        check(f"mini_{bl}_{br}_{tl}_cm", mo == cm_n, detail)
        check(f"mini_{bl}_{br}_{tl}_H_cap_L", hl == 0,
              f"mini ({bl},{br},{tl}): H cap L = {hl}")

    # ------------------------------------------------------------------
    # Chain inequality brute force
    # ------------------------------------------------------------------
    ok_chain = True
    for total in range(2, SERIES_DEPTH + 1):
        for n in range(1, total):
            m = total - n
            if not (c_all[total] <= c_all[m] * (c_all[n] - a_series[n - 1])):
                ok_chain = False
    check("chain_small_instances", ok_chain, "c_(m+n) <= c_m (c_n - a_(n-1)) for m+n<=9")

    # ------------------------------------------------------------------
    # Integer identities, traps, and JSON consistency
    # ------------------------------------------------------------------
    h_union_l = H35 + L35
    union = h_union_l + C35 - c_cap_l - c_cap_h
    M = C36 - union
    check("union_identity", union == 76401062626946721319, f"union = {union}")
    check("M_identity", M == 2941294455272074779839351, f"M = {M}")
    naive = H35 + L35 + C35
    check("trap_no_overlap_subtraction_is_wrong", naive != union and naive - union == c_cap_l + c_cap_h,
          "forgetting the overlap subtraction gives a different (wrong) number")
    check("H_union_L_identity", h_union_l == 76117709690825156721,
          "H cap L = 0, so the exact H/L union has the required audited value")
    check("trap_H_cap_L_zero_is_essential", h_union_l <= union,
          "H cap L = 0 keeps the pairwise H/L sum valid")
    E = C35 - 3 ** 35 - 2 * (35 * 2 ** 34)
    check("E_route_identities",
          E == 283510542706105118
          and E - c_cap_h == 283352901761826230
          and 76117709690825156721 + (E - c_cap_h) == 76401062592586982951
          and C36 - 76401062592586982951 == 2941294455272109139577719,
          "conservative E-route constants (the supplied prompt values, now derived)")
    check("E_route_strictly_weaker",
          76401062592586982951 < union and union - 76401062592586982951 == 4 * 2 ** 33,
          "E-route union is smaller by exactly 4*2^33")
    check("M_beats_wave10", M < C36 - H35 < C36, "M < c_36 - H_35 (wave-10 base)")

    payload = json.loads(RESULT.read_text(encoding="utf-8"))
    data = payload["data"]
    check("json_provenance", payload["provenance"]["script"] == "experiments/e116_saw_union.py"
          and payload["provenance"]["mpmath_iv_dps"] >= 100, "provenance fields")
    check("json_checks_all_passed",
          len(payload["checks"]) > 0 and all(c["passed"] for c in payload["checks"]),
          f"{len(payload['checks'])} producer checks all passed")
    check("json_constants",
          data["family_H"]["B_11_7"] == B
          and data["family_H"]["T_13"] == T
          and data["family_H"]["H_35"] == H35
          and data["family_L"]["L_35"] == L35
          and data["family_C"]["C_35"] == C35
          and data["overlaps"]["C_cap_L"]["value"] == c_cap_l
          and data["overlaps"]["C_cap_H"]["value"] == c_cap_h
          and data["overlaps"]["H_cap_L"]["value"] == 0
          and data["union_certificate"]["components"]["H_union_L"] == h_union_l
          and data["conservative_E_route"]["H_union_L"] == h_union_l
          and data["union_certificate"]["value"] == union
          and data["chain_theorem"]["n36_instance"]["M"] == M,
          "stored JSON constants equal the independently recomputed values")
    check("json_external_c36", data["external_source"]["c_36"] == C36
          and data["external_source"]["sha256"] == C36_SHA256, "stored external c_36 and hash")
    n36 = data["chain_theorem"]["n36_instance"]
    check("json_chain_metadata",
          n36["M_wave10"] == C36 - H35 and n36["M_E_route"] == C36 - 76401062592586982951
          and n36["a_35_lower"] == union, "chain metadata fields consistent")

    # ------------------------------------------------------------------
    # Endpoint enclosure and directed floor
    # ------------------------------------------------------------------
    lo, hi = enclosure(M, 36)
    f = floor40(lo)
    check("endpoint_floor", f == "0.2122120589214465859334619330429416874541",
          f"floor40(atanh(M^(-1/36))) = {f}")
    with localcontext() as ctx:
        ctx.prec = 220
        check("directed_floor_semantics",
              Decimal(f) <= Decimal(lo) and Decimal(hi) < Decimal(f) + Decimal(1).scaleb(-40),
              "floor <= lo <= x <= hi < floor + 10^-40 at 150 dps")
    lo92, _ = enclosure(C36 - H35, 36)
    loE, _ = enclosure(C36 - 76401062592586982951, 36)
    lo0, _ = enclosure(C36, 36)
    check("endpoint_ordering",
          Decimal(floor40(lo0)) < Decimal(floor40(lo92)) < Decimal(floor40(loE)) < Decimal(f),
          "incumbent < wave-10 < E-route < new endpoint floors")
    check("json_endpoint_matches",
          data["endpoint"]["new_floor_40"] == f
          and data["endpoint"]["wave10_floor_40"] == floor40(lo92)
          and data["endpoint"]["E_route_floor_40"] == floor40(loE)
          and data["endpoint"]["incumbent_floor_40"] == floor40(lo0),
          "stored floors equal the 150-dps recomputation")

    print()
    if FAILURES:
        print(f"FAILED ({len(FAILURES)}): {FAILURES}")
        raise SystemExit(1)
    print("OK")


if __name__ == "__main__":
    main()
