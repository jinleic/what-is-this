#!/usr/bin/env python3
"""Standalone verifier for experiments/e144_octahedral.py.

Independently recomputes (no producer import):
  * the order-48 octahedral site action on the 27-site box {0,1,2}^3 and the
    8-site box {0,1}^3, by brute-force orbit counting with canonical forms;
  * the column and anchored-row Burnside counts by an independent DP;
  * the closed forms |C(3,2,s)| and rank S_s;
  * the full-class s<=2 and s<=3 echelon ranks at both primes with a fresh
    implementation of D = (1/2)[H, .] images, translation-orbit anchoring and
    column echelon over F_p;
and checks the artifact results/integrability/octahedral_s6.json against all
of these, plus the decision arithmetic (cap lemma, memory models, verdict).

Run: PYTHONPATH=src .venv/bin/python tests/test_octahedral.py
"""

from __future__ import annotations

import itertools
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "integrability" / "octahedral_s6.json"

PRIMES = (2_147_483_647, 2_147_483_629)

FAILURES: list[str] = []


def check(name: str, passed: bool, detail: str = "") -> None:
    status = "PASS" if passed else "FAIL"
    print(f"[{status}] {name}" + (f": {detail}" if detail else ""))
    if not passed:
        FAILURES.append(name)


# ---------------------------------------------------------------------------
# Independent group and Burnside machinery.
# ---------------------------------------------------------------------------

BASE27 = tuple(itertools.product(range(3), repeat=3))
IDX27 = {s: i for i, s in enumerate(BASE27)}
CTR27 = tuple((a - 1, b - 1, c - 1) for a, b, c in BASE27)


def build_perms_centered(sites):
    """Site permutation tables for the 48 signed permutations about the
    centre of the given box (odd side: integer centre, even side: half)."""
    n = len(sites[0])
    side = max(s[0] for s in sites) + 1
    # doubled centred coordinates: D = 2a - (side-1)
    dbl = lambda s: tuple(2 * s[i] - (side - 1) for i in range(n))
    dmap = {dbl(s): i for i, s in enumerate(sites)}
    perms = []
    for p in itertools.permutations(range(n)):
        for signs in itertools.product((1, -1), repeat=n):
            table = []
            ok = True
            for s in sites:
                d = dbl(s)
                img = tuple(signs[i] * d[p[i]] for i in range(n))
                if img not in dmap:
                    ok = False
                    break
                table.append(dmap[img])
            if ok:
                perms.append(tuple(table))
    return perms


PERMS27 = build_perms_centered(BASE27)


def cycle_lengths(perm):
    seen = [False] * len(perm)
    out = []
    for i in range(len(perm)):
        if seen[i]:
            continue
        j, L = i, 0
        while not seen[j]:
            seen[j] = True
            j = perm[j]
            L += 1
        out.append(L)
    return out


def fix_counts(cycles, smax):
    ways = [1] + [0] * smax
    for L in cycles:
        new = ways[:]
        for n in range(L, smax + 1):
            new[n] += 3 * ways[n - L]
        ways = new
    fix, tot = [], 0
    for n in range(smax + 1):
        tot += ways[n]
        fix.append(tot)
    return fix


def column_orbits(perms, n_sites, smax):
    assert len(set(perms)) == 48, f"expected 48 distinct site actions, got {len(set(perms))}"
    fixes = [fix_counts(cycle_lengths(p), smax) for p in perms]
    orbits = []
    for s in range(smax + 1):
        tot = sum(f[s] for f in fixes)
        check(f"burnside_integrality_s{s}" if n_sites == 27 else f"b1_burnside_integrality_s{s}",
              tot % 48 == 0, f"sum of fixes {tot} divisible by 48")
        orbits.append(tot // 48)
    return orbits


def brute_column_orbits(perms, n_sites, smax):
    canon = set()
    for k in range(smax + 1):
        for sites in itertools.combinations(range(n_sites), k):
            for codes in itertools.product((1, 2, 3), repeat=k):
                w = tuple(zip(sites, codes))
                best = None
                for perm in perms:
                    key = tuple(sorted((perm[i], c) for i, c in w))
                    if best is None or key < best:
                        best = key
                canon.add(best)
    return len(canon)


def affine_row_orbits(smax):
    """Burnside count of G-orbits of anchored words on the 27-site box via
    finite in-box cycles of A_t(y) = g y - t (centred coordinates)."""
    per_g = []
    for p in itertools.permutations(range(3)):
        for signs in itertools.product((1, -1), repeat=3):
            fix = [1] * (smax + 1)  # the identity word
            for t in itertools.product(range(-2, 3), repeat=3):
                cycles, seen_c, done = [], set(), set()
                for start in range(27):
                    if start in done:
                        continue
                    walk, posmap, cur = [start], {start: 0}, start
                    while True:
                        cy = CTR27[cur]
                        ny = tuple(signs[i] * cy[p[i]] - t[i] for i in range(3))
                        if not all(-1 <= v <= 1 for v in ny):
                            break
                        nxt = IDX27[tuple(v + 1 for v in ny)]
                        if nxt in posmap:
                            cyc = frozenset(walk[posmap[nxt]:])
                            if cyc not in seen_c:
                                seen_c.add(cyc)
                                cycles.append(cyc)
                            break
                        posmap[nxt] = len(walk)
                        walk.append(nxt)
                        cur = nxt
                    done.update(walk)
                dp = [[0] * 8 for _ in range(smax + 1)]
                dp[0][0] = 1
                for cyc in cycles:
                    L, cov = len(cyc), 0
                    for i in cyc:
                        c = CTR27[i]
                        if c[0] == -1:
                            cov |= 1
                        if c[1] == -1:
                            cov |= 2
                        if c[2] == -1:
                            cov |= 4
                    new = [row[:] for row in dp]
                    for n in range(L, smax + 1):
                        for c0 in range(8):
                            if dp[n - L][c0]:
                                new[n][c0 | cov] += 3 * dp[n - L][c0]
                    dp = new
                for s in range(smax + 1):
                    fix[s] += sum(dp[n][7] for n in range(1, s + 1))
            per_g.append(fix)
    orbits = []
    for s in range(smax + 1):
        tot = sum(f[s] for f in per_g)
        check(f"row_burnside_integrality_s{s}", tot % 48 == 0)
        orbits.append(tot // 48)
    return orbits


def brute_row_orbits(smax):
    def rekey(pairs):
        pts = [(BASE27[i], c) for i, c in pairs]
        shift = [min(p[0][a] for p in pts) for a in range(3)]
        return tuple(sorted(((p[0][0] - shift[0], p[0][1] - shift[1],
                              p[0][2] - shift[2]), p[1]) for p in pts))

    canon = set()
    for k in range(0, smax + 1):
        for sites in itertools.combinations(range(27), k):
            if k:
                pts = [BASE27[i] for i in sites]
                if not (min(p[0] for p in pts) == 0 and min(p[1] for p in pts) == 0
                        and min(p[2] for p in pts) == 0):
                    continue
            for codes in itertools.product((1, 2, 3), repeat=k):
                w = tuple(zip(sites, codes))
                best = None
                for perm in PERMS27:
                    g = [(perm[i], c) for i, c in w]
                    key = rekey(g) if g else ()
                    if best is None or key < best:
                        best = key
                canon.add(best)
    return len(canon)


# ---------------------------------------------------------------------------
# Independent class-system echelon (small sizes).
# ---------------------------------------------------------------------------

EXT = tuple(itertools.product(range(-1, 4), repeat=3))
EXTI = {s: i for i, s in enumerate(EXT)}
BITS27 = tuple(1 << EXTI[s] for s in BASE27)
NB = {}
for _s, _b in zip(BASE27, BITS27):
    _lst = []
    for _ax in range(3):
        for _sg in (-1, 1):
            _v = list(_s)
            _v[_ax] += _sg
            _lst.append(1 << EXTI[tuple(_v)])
    NB[_b] = tuple(_lst)


def image(x, z):
    out = {}
    m = z
    while m:
        bit = m & -m
        m ^= bit
        k = (x ^ bit, z)
        out[k] = out.get(k, 0) + 1
    m = x
    while m:
        bit = m & -m
        m ^= bit
        for nb in NB[bit]:
            if not x & nb:
                k = (x, z ^ bit ^ nb)
                out[k] = out.get(k, 0) - 1
    return {k: v for k, v in out.items() if v}


def okey(x, z):
    support = x | z
    if not support:
        return 0
    mnx = mny = mnz = 99
    pts = []
    m, size = support, 0
    while m:
        bit = m & -m
        m ^= bit
        pos = bit.bit_length() - 1
        sx, sy, sz = EXT[pos]
        mnx, mny, mnz = min(mnx, sx), min(mny, sy), min(mnz, sz)
        pts.append((pos, sx, sy, sz))
        size += 1
    digs = []
    for pos, sx, sy, sz in pts:
        code = (1 if (x >> pos) & 1 else 0) | (2 if (z >> pos) & 1 else 0)
        digs.append((((sx - mnx) * 5 + (sy - mny)) * 5 + (sz - mnz)) * 3 + (code - 1))
    digs.sort()
    packed = 0
    for d in digs:
        packed = packed * 375 + d
    return (size << 60) | packed


def class_columns(smax):
    yield (0, 0)
    for size in range(1, smax + 1):
        for pos in itertools.combinations(range(27), size):
            bits = [BITS27[i] for i in pos]
            for codes in itertools.product((1, 2, 3), repeat=size):
                x = z = 0
                for b, code in zip(bits, codes):
                    if code & 1:
                        x |= b
                    if code & 2:
                        z |= b
                yield (x, z)


def echelon_rank(smax, prime):
    pivots: dict[int, dict[int, int]] = {}
    for x, z in class_columns(smax):
        w = {}
        for (nx, nz), c in image(x, z).items():
            k = okey(nx, nz)
            w[k] = (w.get(k, 0) + c) % prime
        w = {k: v for k, v in w.items() if v}
        while w:
            lead = max(w)
            prior = pivots.get(lead)
            if prior is None:
                inv = pow(w[lead], prime - 2, prime)
                pivots[lead] = {k: v * inv % prime for k, v in w.items() if k != lead}
                break
            sc = w.pop(lead)
            for k, v in prior.items():
                nv = (w.get(k, 0) - sc * v) % prime
                if nv:
                    w[k] = nv
                else:
                    w.pop(k, None)
    return len(pivots)


# ---------------------------------------------------------------------------
# Main.
# ---------------------------------------------------------------------------

def main() -> int:
    print("test_octahedral: independent verification of e144_octahedral\n")

    # 1. Group + closed forms.
    check("group_has_48_distinct_site_actions", len(set(PERMS27)) == 48)
    cols_by_s = [sum(math.comb(27, k) * 3 ** k for k in range(s + 1)) for s in range(7)]
    check("columns_s6_closed_form", cols_by_s[6] == 236_912_446, f"{cols_by_s[6]:,}")
    check("columns_s5_matches_frozen", cols_by_s[5] == 21_121_156)
    A = [math.comb(27, k) - 3 * math.comb(18, k) + 3 * math.comb(12, k) - math.comb(8, k)
         for k in range(1, 7)]
    rankS = [1 + sum(A[k - 1] * 3 ** k for k in range(1, s + 1)) for s in range(7)]
    check("rankS_s4_matches_frozen", rankS[4] == 822_334)
    check("rankS_s5_matches_frozen", rankS[5] == 14_757_412)
    check("rankS_s6_value", rankS[6] == 191_949_610, f"{rankS[6]:,}")

    # 2. Burnside vs brute force.
    print("\n[brute-force orbit counting]")
    orbits27 = column_orbits(PERMS27, 27, 6)
    bf2 = brute_column_orbits(PERMS27, 27, 2)
    check("column_burnside_vs_bruteforce_s2", bf2 == orbits27[2],
          f"brute {bf2} vs DP {orbits27[2]}")
    rorbits = affine_row_orbits(6)
    bfrows2 = brute_row_orbits(2)
    check("row_burnside_vs_bruteforce_s2", bfrows2 == rorbits[2],
          f"brute {bfrows2} vs DP {rorbits[2]}")

    B1 = tuple(itertools.product(range(2), repeat=3))
    PERMS8 = build_perms_centered(B1)
    orbits8 = column_orbits(PERMS8, 8, 8)
    bf8 = brute_column_orbits(PERMS8, 8, 8)
    check("b1_column_burnside_vs_bruteforce_s8", bf8 == orbits8[8],
          f"brute {bf8} vs DP {orbits8[8]}")

    # 3. Independent echelons.
    print("\n[independent small-class echelons]")
    r2 = [echelon_rank(2, p) for p in PRIMES]
    check("echelon_s2_both_primes", r2[0] == r2[1] == rankS[2] - 2,
          f"ranks {r2} vs rank S - 2 = {rankS[2] - 2}")
    r3 = [echelon_rank(3, p) for p in PRIMES]
    check("echelon_s3_both_primes", r3[0] == r3[1] == rankS[3] - 2,
          f"ranks {r3} vs rank S - 2 = {rankS[3] - 2}")

    # 4. Artifact consistency.
    print("\n[artifact checks]")
    if not ARTIFACT.exists():
        check("artifact_exists", False, str(ARTIFACT))
        print(f"\n{len(FAILURES)} failing checks")
        return 1
    art = json.loads(ARTIFACT.read_text())
    data = art["data"]
    feas = data["feasibility"]

    check("artifact_all_checks_passed",
          all(c["passed"] for c in art["checks"]),
          f"{sum(c['passed'] for c in art['checks'])}/{len(art['checks'])}")

    cb = feas["column_burnside"]
    check("artifact_column_orbits_match_independent_dp",
          cb["orbit_counts_by_s"] == orbits27,
          f"artifact {cb['orbit_counts_by_s'][6]:,} vs recomputed {orbits27[6]:,}")
    check("artifact_column_orbits_s6_literal",
          cb["orbit_counts_by_s"][6] == 5_048_368)
    rb = feas["row_burnside"]
    check("artifact_row_orbits_match_independent_dp",
          rb["orbit_counts_by_s"] == rorbits)
    check("artifact_row_orbits_s6_literal",
          rb["orbit_counts_by_s"][6] == 4_094_477)
    check("artifact_fix_sums_divisible_by_48",
          all(v % 48 == 0 for v in cb["per_element_fix_sum_by_s"])
          and all(v % 48 == 0 for v in rb["per_element_fix_sum_by_s"]))
    check("artifact_columns_s6", feas["columns_s6"] == cols_by_s[6] == 236_912_446)
    check("artifact_rankS_by_s", feas["rank_S_by_s"] == rankS)
    check("artifact_b1_column_orbits_s8", feas["b1_box_machine_check"]["column_orbits_by_s"][8] == bf8)
    check("artifact_b1_row_orbits_s4", feas["b1_box_machine_check"]["row_orbits_by_s"][4] == 269)

    reg = data["regression_full_class"]
    s2r = [r.get("rank_Fp") for r in reg["s2"]["runs"] if not r.get("walled")]
    s3r = [r.get("rank_Fp") for r in reg["s3"]["runs"] if not r.get("walled")]
    s4r = [r.get("rank_Fp") for r in reg["s4_frozen_anchor"]["runs"] if not r.get("walled")]
    check("artifact_s2_ranks_match_independent_echelon",
          s2r == r2, f"artifact {s2r} vs independent {r2}")
    check("artifact_s3_ranks_match_independent_echelon",
          s3r == r3, f"artifact {s3r} vs independent {r3}")
    check("artifact_s4_regression_frozen", s4r == [822_332] * 2, f"{s4r}")

    rep5 = data["rep_calibration_s5"]
    rep6 = data["headline_s6_reps"]
    check("artifact_rep5_columns_are_orbits", rep5["columns"] == orbits27[5] == 464_957)
    check("artifact_rep6_columns_are_orbits", rep6["columns"] == orbits27[6] == 5_048_368)
    r5 = [r.get("rank_Fp") for r in rep5["runs"] if not r.get("walled")]
    r6 = [r.get("rank_Fp") for r in rep6["runs"] if not r.get("walled")]
    check("artifact_rep5_ranks_within_cap",
          len(r5) == 2 and all(r <= orbits27[5] for r in r5), f"{r5}")
    check("artifact_rep6_ranks_within_cap",
          len(r6) == 2 and all(r <= orbits27[6] for r in r6), f"{r6}")
    check("artifact_certified_bound_is_min_prime_rank",
          rep6["certified_lower_bound_rank_Q_M6"] == min(r6), f"{r6}")
    d5 = {r.get("pivot_trace_sha256_of_column_lead_stream") for r in rep5["runs"]
          if not r.get("walled")}
    d6 = {r.get("pivot_trace_sha256_of_column_lead_stream") for r in rep6["runs"]
          if not r.get("walled")}
    print(f"  [obs] rep lead-stream digest prime agreement: s5={'yes' if len(d5) == 1 else 'no'}, "
          f"s6={'yes' if len(d6) == 1 else 'no'}")

    dec = data["decision"]
    check("artifact_required_pivots", dec["required_pivots"] == rankS[6] - 2 == 191_949_608)
    check("artifact_cap_lemma_holds",
          dec["cap_lemma"]["orbit_representatives"] < dec["required_pivots"]
          and dec["cap_lemma"]["orbit_representatives"] == orbits27[6])
    mm = dec["memory_models_bytes"]
    check("artifact_memory_models_exceed_8gib",
          mm["total_lead_plus_csr"] > 8 * 1024 ** 3
          and mm["whole_process_model"] > 8 * 1024 ** 3,
          f"lead+CSR {mm['total_lead_plus_csr']:,} B, whole-process {mm['whole_process_model']:,} B")
    check("artifact_lead_plus_csr_model_recomputes",
          mm["csr_pivot_store"] == (191_949_608 * 313_579_672 * 12) // 14_757_410
          and mm["lead_table"] == 12 << 28)
    check("artifact_did_not_decide_s6",
          rep6["certified_lower_bound_rank_Q_M6"] < dec["required_pivots"],
          "representative lower bound strictly below required pivots")

    print(f"\n{len(FAILURES)} failing checks" + ("" if not FAILURES else f": {FAILURES}"))
    return 1 if FAILURES else 0


if __name__ == "__main__":
    raise SystemExit(main())
