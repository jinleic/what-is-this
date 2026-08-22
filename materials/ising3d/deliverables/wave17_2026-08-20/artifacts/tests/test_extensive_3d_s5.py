#!/usr/bin/env python3
"""Independent verifier for the Z^3 R=2 support-class s<=5 census (e135).

This test deliberately does not import `experiments/e135_extensive_3d_s5.py`.
Everything decisive is rebuilt from raw inputs with separate representations:

* the producer's orbit rows are single integers (size<<54 | base-500 digits);
  this verifier uses a distinct bytes encoding (sorted relative-site bytes +
  one ordered-Pauli code byte per site) and a plain dict lead map;
* the producer certifies at primes 2^31-1 and 2^31-19; this verifier
  eliminates at the independently proved third prime 1,000,000,007;
* the generator convention itself (D = (1/2)[H, .], ordered X^a Z^b, a=b=1)
  is cross-validated against exact dense matrix commutators on a seven-site
  closed neighbourhood, not only against the producer's ad formulas.

Rebuilt independently:
  1. dense-matrix commutator spot check of the ad convention;
  2. the closed-form ranking arithmetic: class column counts, rank S by
     brute anchored-subset enumeration AND by inclusion-exclusion, the
     analytic target rank S - 2, and the kernel/nullity arithmetic;
  3. the full s<=4 regression rank 822,332 (plus the anchored census
     822,334 and the s<=3 anchors 82,216 / 29,749 / 29,747) at the third
     prime with the independent encoding;
  4. the producer JSON consistency gate: a quotient may be claimed ONLY
     for cases with both primes completed inside their predeclared
     process_time budgets and the RSS wall; any wall on a claimed case,
     any budget breach, or any arithmetic mismatch FAILS this test.  An
     honestly recorded OBSERVED wall for the s<=5 class passes the gate
     only when no certified claim is attached to it.

Run:
    PYTHONPATH=src .venv/bin/python tests/test_extensive_3d_s5.py
"""
from __future__ import annotations

import itertools
import json
import math
import platform
import resource
import time
from pathlib import Path
from typing import Iterator

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "integrability" / "extensive_3d_s5.json"

PRODUCER_PRIMES = (2_147_483_647, 2_147_483_629)
TEST_PRIME = 1_000_000_007
# Predeclared verifier budgets (process_time; a breach is a FAILURE, since
# every value rebuilt here is claimed).
S4_ELIMINATION_BUDGET_S = 1_500.0
S3_ELIMINATION_BUDGET_S = 120.0
ARITHMETIC_BUDGET_S = 600.0
RSS_WALL_BYTES = 6_000_000_000  # the producer's predeclared wall

CERTIFIED = {
    "s2_columns": 3_241, "s2_rank_S": 562, "s2_rank": 560,
    "s3_columns": 82_216, "s3_rank_S": 29_749, "s3_rank": 29_747,
    "s4_columns": 1_503_766, "s4_rank_S": 822_334, "s4_rank": 822_332,
    "s5_columns": 21_121_156, "s5_rank_S": 14_757_412, "s5_rank_target": 14_757_410,
}

FAILURES: list[str] = []


def require(name: str, condition: bool, detail: str = "") -> None:
    status = "ok" if condition else "FAIL"
    print(f"  [{status}] {name}" + (f" — {detail}" if detail else ""), flush=True)
    if not condition:
        FAILURES.append(name)


def peak_rss_bytes() -> int:
    raw = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return raw if platform.system() == "Darwin" else raw * 1024

def is_prime_by_trial_division(value: int) -> bool:
    if value < 2:
        return False
    if value % 2 == 0:
        return value == 2
    d = 3
    while d * d <= value:
        if value % d == 0:
            return False
        d += 2
    return True


# ---------------------------------------------------------------------------
# Independent model: bytes orbit encoding.
# ---------------------------------------------------------------------------

BASE_SITES = tuple(itertools.product(range(3), repeat=3))
EXT_SITES = tuple(itertools.product(range(-1, 4), repeat=3))
EXT_INDEX = {site: i for i, site in enumerate(EXT_SITES)}
BASE_BITS = tuple(1 << EXT_INDEX[s] for s in BASE_SITES)
EXT_MASK = (1 << len(EXT_SITES)) - 1
NEIGH: dict[int, tuple[int, ...]] = {}
for _s, _b in zip(BASE_SITES, BASE_BITS):
    _ns = []
    for _ax in range(3):
        for _sg in (-1, 1):
            _n = list(_s)
            _n[_ax] += _sg
            _ns.append(1 << EXT_INDEX[tuple(_n)])
    NEIGH[_b] = tuple(_ns)


def image_terms(x: int, z: int) -> dict[tuple[int, int], int]:
    """(1/2)[H, X^x Z^z] by the ad formulas (validated against matrices)."""
    out: dict[tuple[int, int], int] = {}
    m = z
    while m:
        b = m & -m
        m ^= b
        k = (x ^ b, z)
        out[k] = out.get(k, 0) + 1
    m = x
    while m:
        b = m & -m
        m ^= b
        for nb in NEIGH[b]:
            if not x & nb:
                k = (x, z ^ b ^ nb)
                out[k] = out.get(k, 0) - 1
    return {k: v for k, v in out.items() if v}


def canonical(x: int, z: int) -> bytes:
    """Distinct bytes key per translation orbit: sorted support after
    subtracting the coordinatewise minimum, 1 byte per relative coordinate
    and 1 byte for the local ordered-Pauli code (X=1, Z=2, XZ=3)."""
    support = x | z
    if not support:
        return b""
    sites = []
    m = support
    while m:
        b = m & -m
        m ^= b
        sites.append(EXT_SITES[b.bit_length() - 1])
    mn = [min(s[a] for s in sites) for a in range(3)]
    entries = []
    for site in sites:
        rel = tuple(site[a] - mn[a] for a in range(3))
        bit = 1 << EXT_INDEX[site]
        code = (1 if x & bit else 0) | (2 if z & bit else 0)
        entries.append(rel + (code,))
    entries.sort()
    return b"".join(bytes(e) for e in entries)


def explicit_kernel_traps() -> dict[str, bool]:
    """Independent orbit-sum checks for D(Z_0), D(X_0), and pi D h."""
    origin = 1 << EXT_INDEX[(0, 0, 0)]
    dz_ok = image_terms(0, origin) == {(origin, origin): 1}
    dx = image_terms(origin, 0)
    dx_ok = len(dx) == 6 and set(dx.values()) == {-1}
    h_words: dict[tuple[int, int], int] = {(origin, 0): 1}
    for axis in range(3):
        s = [0, 0, 0]
        s[axis] = 1
        nb = 1 << EXT_INDEX[tuple(s)]
        h_words[(0, origin | nb)] = 1
    orbit_sum: dict[bytes, int] = {}
    for (x, z), coefficient in h_words.items():
        for (nx, nz), c in image_terms(x, z).items():
            k = canonical(nx, nz)
            orbit_sum[k] = orbit_sum.get(k, 0) + coefficient * c
    return {
        "D_Z0": dz_ok,
        "D_X0": dx_ok,
        "pi_D_h": not any(orbit_sum.values()),
    }

def iter_class_columns(smax: int) -> Iterator[tuple[int, int]]:
    yield (0, 0)
    for size in range(1, smax + 1):
        for positions in itertools.combinations(range(27), size):
            bits = [BASE_BITS[p] for p in positions]
            for codes in itertools.product((1, 2, 3), repeat=size):
                x = z = 0
                for b, code in zip(bits, codes):
                    if code & 1:
                        x |= b
                    if code & 2:
                        z |= b
                yield (x, z)


def is_anchored(support: int) -> bool:
    if not support:
        return True
    mn = [9] * 3
    m = support
    while m:
        b = m & -m
        m ^= b
        s = EXT_SITES[b.bit_length() - 1]
        for a in range(3):
            if s[a] < mn[a]:
                mn[a] = s[a]
    return mn == [0, 0, 0]


def eliminate_bytes(smax: int, prime: int, budget_s: float) -> dict[str, int | float]:
    """Independent column echelon over F_p with bytes orbit keys.

    Any total order on keys gives the same pivot count (the rank), so the
    plain bytes order is used.  Returns counters for the assertions.
    """
    t0 = time.process_time()
    pivots: dict[bytes, list[tuple[bytes, int]]] = {}
    rank = reductions = anchored = 0
    cols = 0
    for x, z in iter_class_columns(smax):
        cols += 1
        if is_anchored(x | z):
            anchored += 1
        w: dict[bytes, int] = {}
        for (nx, nz), c in image_terms(x, z).items():
            assert not (nx | nz) & ~EXT_MASK
            assert (nx | nz).bit_count() <= smax + 1
            k = canonical(nx, nz)
            w[k] = (w.get(k, 0) + c) % prime
        w = {k: v for k, v in w.items() if v}
        while w:
            lead = max(w)
            pv = pivots.get(lead)
            if pv is None:
                inv = pow(w[lead], prime - 2, prime)
                vec = sorted((k, v * inv % prime) for k, v in w.items() if k != lead)
                pivots[lead] = vec
                rank += 1
                break
            reductions += 1
            sc = w.pop(lead)
            for k, v in pv:
                nv = (w.get(k, 0) - sc * v) % prime
                if nv:
                    w[k] = nv
                else:
                    w.pop(k, None)
        if cols % 500_000 == 0:
            elapsed = time.process_time() - t0
            rss = peak_rss_bytes()
            if elapsed > budget_s or rss > RSS_WALL_BYTES:
                raise RuntimeError(
                    f"verifier wall at {cols} columns: process_time={elapsed:.3f}s "
                    f"(budget {budget_s}s), peak RSS={rss} bytes "
                    f"(wall {RSS_WALL_BYTES})")
    elapsed = time.process_time() - t0
    rss = peak_rss_bytes()
    if elapsed > budget_s or rss > RSS_WALL_BYTES:
        raise RuntimeError(
            f"verifier wall at stage end: process_time={elapsed:.3f}s "
            f"(budget {budget_s}s), peak RSS={rss} bytes "
            f"(wall {RSS_WALL_BYTES})")
    return {"columns": cols, "rank": rank, "reductions": reductions,
            "anchored": anchored, "elapsed": round(elapsed, 3), "peak_rss_bytes": rss}


# ---------------------------------------------------------------------------
# Dense-matrix cross-validation of the generator convention over Q.
# ---------------------------------------------------------------------------

def matrix_commutor_check() -> bool:
    """Window-internal validation of the ad convention against dense exact
    matrix commutators over Z (int64 numpy; all entries are small integers,
    so the arithmetic is exact with no overflow anywhere near 2^63).

    Window W = {(0,0,1)} plus its six neighbours (7 sites, dim 128).
    H_W = sum_{s in W} X_s + sum_{window-internal bonds uv} Z_u Z_v.
    For each tested ordered-Pauli word w supported inside W,

        (1/2)[H_W, w]  ==  sum of ad-image terms of w with support in W

    as an exact matrix identity: ad-image terms leaving W arise only from
    bonds crossing the window boundary (their output words contain the
    outside endpoint), and operators disjoint from w commute, so the
    window-internal part of the full commutator is exactly [H_W, w].
    Tested for all 255 nonidentity words on the 4-site subwindow
    {(0,0,1),(1,0,1),(0,1,1),(0,0,2)} (all in B_2, all bonded to the
    centre) and all 21 single-site words of W.
    """
    import numpy as np

    centre = (0, 0, 1)
    window = [centre] + [tuple(centre[a] + d[a] for a in range(3))
                         for d in ((1, 0, 0), (-1, 0, 0), (0, 1, 0),
                                   (0, -1, 0), (0, 0, 1), (0, 0, -1))]
    wbits = [1 << EXT_INDEX[s] for s in window]
    covered = 0
    for b in wbits:
        covered |= b
    n = len(window)
    dim = 1 << n

    x1 = np.array([[0, 1], [1, 0]], dtype=np.int64)
    z1 = np.array([[1, 0], [0, -1]], dtype=np.int64)
    eye = np.eye(2, dtype=np.int64)

    def site_op(single, site):
        out = None
        for s in window:
            m = single if s == site else eye
            out = m if out is None else np.kron(out, m)
        return out

    x_ops = {s: site_op(x1, s) for s in window}
    z_ops = {s: site_op(z1, s) for s in window}
    h_w = np.zeros((dim, dim), dtype=np.int64)
    for s in window:
        h_w += x_ops[s]
    for i in range(n):
        for j in range(i + 1, n):
            a, b = window[i], window[j]
            if sum(abs(a[k] - b[k]) for k in range(3)) == 1:
                h_w += z_ops[a] @ z_ops[b]

    def word_matrix(x: int, z: int):
        m = np.eye(dim, dtype=np.int64)
        for i in range(n):
            bit = wbits[i]
            code = (1 if x & bit else 0) | (2 if z & bit else 0)
            if code == 1:
                m = m @ x_ops[window[i]]
            elif code == 2:
                m = m @ z_ops[window[i]]
            elif code == 3:
                m = m @ (x_ops[window[i]] @ z_ops[window[i]])
        return m

    subwindow = [(0, 0, 1), (1, 0, 1), (0, 1, 1), (0, 0, 2)]
    sub_bits = [1 << EXT_INDEX[s] for s in subwindow]
    test_words: list[tuple[int, int]] = []
    for size in range(1, len(sub_bits) + 1):
        for positions in itertools.combinations(range(len(sub_bits)), size):
            for codes in itertools.product((1, 2, 3), repeat=size):
                x = z = 0
                for pos, code in zip(positions, codes):
                    if code & 1:
                        x |= sub_bits[pos]
                    if code & 2:
                        z |= sub_bits[pos]
                test_words.append((x, z))
    for i in range(n):  # every valid base-box single-site word of W
        if window[i] not in BASE_SITES:
            continue
        for code in (1, 2, 3):
            x = wbits[i] if code & 1 else 0
            z = wbits[i] if code & 2 else 0
            test_words.append((x, z))
    for x, z in test_words:
        w_mat = word_matrix(x, z)
        comm = h_w @ w_mat - w_mat @ h_w
        predicted = np.zeros((dim, dim), dtype=np.int64)
        for (nx, nz), c in image_terms(x, z).items():
            assert not (nx | nz) & ~EXT_MASK
            if (nx | nz) & ~covered:
                continue
            predicted += c * word_matrix(nx, nz)
        if not np.array_equal(comm, 2 * predicted):
            print(f"    matrix mismatch at word ({x:#x},{z:#x})")
            return False
    return True


# ---------------------------------------------------------------------------
# Ranking arithmetic, rebuilt.
# ---------------------------------------------------------------------------

def anchored_subset_counts_brute(smax: int) -> list[int]:
    counts = []
    for k in range(1, smax + 1):
        n = 0
        for combo in itertools.combinations(range(27), k):
            sites = [BASE_SITES[c] for c in combo]
            if (min(s[0] for s in sites) == 0 and min(s[1] for s in sites) == 0
                    and min(s[2] for s in sites) == 0):
                n += 1
        counts.append(n)
    return counts


def anchored_subset_count_ie(k: int) -> int:
    return (math.comb(27, k) - 3 * math.comb(18, k)
            + 3 * math.comb(12, k) - math.comb(8, k))


def rank_S(smax: int) -> int:
    brute = anchored_subset_counts_brute(smax)
    assert all(anchored_subset_count_ie(k + 1) == brute[k] for k in range(smax))
    return 1 + sum(a * 3 ** (k + 1) for k, a in enumerate(brute))


def column_count(smax: int) -> int:
    return sum(math.comb(27, k) * 3**k for k in range(smax + 1))


def main() -> int:
    print("== independent verifier: extensive_3d_s5 ==", flush=True)
    t_all = time.process_time()

    require("producer_primes_are_prime",
            all(is_prime_by_trial_division(p) for p in PRODUCER_PRIMES),
            "2,147,483,647 and 2,147,483,629 by deterministic trial division")
    require("third_prime_is_prime", is_prime_by_trial_division(TEST_PRIME),
            f"{TEST_PRIME} by deterministic trial division")

    # -- 1. dense-matrix commutator cross-validation -----------------------
    require("ad_convention_matches_dense_matrix_commutators",
            matrix_commutor_check(),
            "(1/2)[H_W, w] equals the window-internal ad image as an exact "
            "matrix identity on a seven-site closed neighbourhood")

    # -- 2. ranking arithmetic --------------------------------------------
    t0 = time.process_time()
    require("column_counts_closed_form",
            column_count(3) == CERTIFIED["s3_columns"]
            and column_count(4) == CERTIFIED["s4_columns"]
            and column_count(5) == CERTIFIED["s5_columns"],
            "82,216 / 1,503,766 / 21,121,156 = sum C(27,k) 3^k")
    require("rank_S_bruteforce_equals_inclusion_exclusion",
            rank_S(3) == CERTIFIED["s3_rank_S"] and rank_S(4) == CERTIFIED["s4_rank_S"]
            and rank_S(5) == CERTIFIED["s5_rank_S"],
            "29,749 / 822,334 / 14,757,412 = 1 + sum A_k 3^k")
    require("analytic_target_is_rank_S_minus_2",
            CERTIFIED["s5_rank_S"] - 2 == CERTIFIED["s5_rank_target"]
            and CERTIFIED["s4_rank_S"] - 2 == CERTIFIED["s4_rank"],
            "the {I,h} kernel upper bound on the commutator rank")
    require("arithmetic_within_budget", time.process_time() - t0 <= ARITHMETIC_BUDGET_S)

    # -- 3. independent eliminations at the third prime --------------------
    r2 = eliminate_bytes(2, TEST_PRIME, 60.0)
    traps = explicit_kernel_traps()
    require("explicit_I_h_kernel_and_sign_traps", all(traps.values()),
            "D(Z_0)=+X_0Z_0, D(X_0) has six -1 bond terms, and pi D h=0 "
            "under the independent bytes orbit encoding")
    require("s2_rank_rebuilt", r2["rank"] == CERTIFIED["s2_rank"]
            and r2["columns"] == CERTIFIED["s2_columns"]
            and r2["anchored"] == CERTIFIED["s2_rank_S"],
            f"rank {r2['rank']}, columns {r2['columns']}, anchored {r2['anchored']}")
    r3 = eliminate_bytes(3, TEST_PRIME, S3_ELIMINATION_BUDGET_S)
    require("s3_rank_rebuilt", r3["rank"] == CERTIFIED["s3_rank"]
            and r3["columns"] == CERTIFIED["s3_columns"]
            and r3["anchored"] == CERTIFIED["s3_rank_S"],
            f"rank {r3['rank']} in {r3['elapsed']}s, anchored {r3['anchored']}")

    r4 = eliminate_bytes(4, TEST_PRIME, S4_ELIMINATION_BUDGET_S)
    require("s4_regression_rank_rebuilt", r4["rank"] == CERTIFIED["s4_rank"],
            f"822,332 reproduced independently at p={TEST_PRIME} "
            f"in {r4['elapsed']}s process time")
    require("s4_regression_columns_and_anchored_census",
            r4["columns"] == CERTIFIED["s4_columns"]
            and r4["anchored"] == CERTIFIED["s4_rank_S"],
            "1,503,766 columns, anchored census 822,334 over the full stream")
    require("s4_kernel_nullity_arithmetic",
            r4["columns"] - r4["rank"] == 681_434
            and r4["columns"] - CERTIFIED["s4_rank_S"] == 681_432,
            "kernel nullity 681,434 = trivial divergence 681,432 + span{I,h}")

    # -- 4. producer JSON consistency gate ---------------------------------
    if not RESULT.exists():
        require("producer_result_json_present", False, str(RESULT))
    else:
        payload = json.loads(RESULT.read_text(encoding="utf-8"))
        data = payload.get("data", {})
        headline = data.get("headline_s5", {})
        outcome = headline.get("outcome")
        require("producer_outcome_known",
                outcome in {"certified_quotient_1", "observed_wall",
                            "modular_rank_below_analytic_bound", "incomplete",
                            "prime_disagreement"}, f"outcome={outcome}")
        cases = headline.get("cases", [])
        wall = headline.get("observed_wall")

        claimed = outcome == "certified_quotient_1"
        if claimed:
            require("claimed_only_with_both_primes_at_target",
                    len(cases) == 2
                    and all(c["rank_Fp"] == CERTIFIED["s5_rank_target"] for c in cases)
                    and headline.get("certified_quotient_dimension") == 1,
                    "both primes completed at rank S - 2 = 14,757,410")
            require("claimed_case_columns_and_census",
                    all(c["columns"] == CERTIFIED["s5_columns"]
                        and c["anchored_columns"] == CERTIFIED["s5_rank_S"]
                        for c in cases),
                    "21,121,156 columns, anchored census 14,757,412")
            require("claimed_case_kernel_arithmetic",
                    headline.get("kernel_nullity") == CERTIFIED["s5_columns"] - CERTIFIED["s5_rank_target"]
                    and headline.get("trivial_divergence_dimension")
                    == CERTIFIED["s5_columns"] - CERTIFIED["s5_rank_S"],
                    "nullity 6,363,746 = trivial 6,363,744 + span{I,h}")
        # any wall on a claimed case, or any completed stage outside its
        # budget or the RSS wall, FAILS unconditionally
        stages_ok = all(
            c["elapsed_seconds_process_time"] <= c["budget_seconds_process_time"]
            and c["peak_rss_bytes"] <= RSS_WALL_BYTES
            for c in cases)
        require("no_budget_or_rss_breach_on_completed_stages", stages_ok,
                "every completed s<=5 stage inside its predeclared process_time "
                "budget and the 6 GB RSS wall")
        require("no_wall_on_claimed_case", not (claimed and wall is not None),
                "a certified claim carries no wall record")
        require("no_certificate_attached_to_observed_wall",
                wall is None or (not claimed and "certified_quotient_dimension" not in headline),
                "an observed wall has no rank or quotient certificate attached")
        if wall is not None:
            require("wall_record_arithmetic_sane",
                    0 <= wall["processed_columns"] <= wall["total_columns"]
                    == CERTIFIED["s5_columns"] and wall["elapsed_seconds_process_time"] > 0,
                    "processed <= total = 21,121,156, positive process_time")

        reg = data.get("regression_s4", {})
        require("producer_regression_ran_before_headline",
                reg.get("reproduced") is True
                and [c["rank_Fp"] for c in reg.get("cases", [])]
                == [CERTIFIED["s4_rank"], CERTIFIED["s4_rank"]],
                "the producer reproduced 822,332 at both primes before the s<=5 run")
        reg_cases = reg.get("cases", [])
        require("no_wall_or_budget_breach_on_claimed_s4_regression",
                reg.get("observed_wall") is None
                and len(reg_cases) == 2
                and all(c["elapsed_seconds_process_time"] <= c["budget_seconds_process_time"]
                        and c["peak_rss_bytes"] <= RSS_WALL_BYTES for c in reg_cases),
                "both claimed s<=4 regression stages completed inside their process-time "
                "budgets and the 6 GB RSS wall")
        reported_checks = payload.get("checks", [])
        require("producer_reported_checks_all_pass",
                bool(reported_checks) and all(c.get("passed") is True for c in reported_checks),
                "every producer self-check recorded as passed")

        greedy = data.get("greedy_gap", {})
        require("producer_greedy_gap_recorded",
                greedy.get("best_greedy_lower_bound") is not None
                and greedy["best_greedy_lower_bound"] < CERTIFIED["s4_rank"]
                and greedy.get("best_greedy_ratio_of_exact", 1.0) < 1.0,
                "greedy static bound stalls strictly below the exact rank")

        unresolved_text = json.dumps(data.get("unresolved", []))
        require("full_box_z3_r2_stated_undecided",
                "FULL Z^3 radius-2" in unresolved_text and "NOT decided" in unresolved_text,
                "the JSON restates that the full 3x3x3 box is not decided")

    total = time.process_time() - t_all
    print(f"verifier process_time: {total:.1f}s", flush=True)
    if FAILURES:
        print(f"FAIL ({len(FAILURES)}): " + ", ".join(FAILURES))
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
