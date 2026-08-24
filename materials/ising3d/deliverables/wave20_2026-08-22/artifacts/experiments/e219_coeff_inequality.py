r"""e219 -- Concatenation-inequality audit for the exact susceptibility coefficients.

Inputs: the theorem-grade coefficients a_0..a_N derived from finite graphs by
experiments/e218_susceptibility_coeffs.py (read from the shared artifact, not
imported), plus the clearly tagged [EXTERNAL] published values through v^32.

Audited statements (proofs/susceptibility_coeffs.md):

  (S1) submultiplicativity      a_{m+n} <= a_m a_n        (m, n >= 1)
  (S2) supermultiplicativity    a_{m+n} >= a_m a_n        (m, n >= 1)
  (S3) log-concavity            a_n^2 >= a_{n-1} a_{n+1}  (n >= 1)
  (S4) ratio ceiling            a_{n+1} <= 6 a_n
  (S5) eventual supermultiplicativity anchored at k: a_{2k} >= a_k^2

Findings recorded here as computed checks:
  * (S2) fails at the smallest possible pair (m, n) = (1, 1):
    a_2 = 30 < 36 = a_1^2, theorem-grade (both values are re-derived and
    hand-countable).  Every anchored variant (S5) fails at every testable k.
  * (S1), (S3), (S4) hold at every checkable index; all remain finite
    evidence only and are tagged as conjectures in the proof note.
  * The natural coefficient-monotonicity / nonnegativity bridge premises are
    FALSE: on the 2x2x1 corner box, [v^5] chi_box = -2 < 0, and the
    coefficient decreases when the box is enlarged from 2x1x1 to 2x2x1.
    The absolute-value domination |[v^n] chi_box| <= a_n survives on every
    tested box (needed by the conditional bridge lemma in the proof note).
  * The plaquette-dressing obstruction: the order-4 dressing of a
    concatenated pair of edge clusters factorizes exactly at distance > 4
    and fails for adjacent clusters (-28 vs -40), which is precisely why the
    SAW concatenation proof does not transfer to dressed coefficients.

Box coefficients here are computed by a vertex-retirement transfer DP over
edge subsets -- an algorithm disjoint from e218's cluster expansion -- and
validated against complete brute-force enumeration on two boxes.

Exact integer arithmetic throughout; single process; trivial memory.
"""

from __future__ import annotations

import json
import resource
import time
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path

SCRIPT = "experiments/e219_coeff_inequality.py"
ROOT = Path(__file__).resolve().parents[1]
RESULT_PATH = ROOT / "results" / "bounds" / "susceptibility_coeffs.json"

K_BOX = 8  # box-study truncation order

# ----------------------------------------------------------------------
# packed geometry (same conventions as e218, re-stated independently)
# ----------------------------------------------------------------------
OFF = 12
STEP = (1 << 10, 1 << 5, 1)


def pack(x: int, y: int, z: int) -> int:
    return ((x + OFF) << 10) | ((y + OFF) << 5) | (z + OFF)


def box_vertices(dims):
    return [
        pack(x, y, z)
        for x in range(dims[0])
        for y in range(dims[1])
        for z in range(dims[2])
    ]


def box_edges(dims):
    es = []
    for x in range(dims[0]):
        for y in range(dims[1]):
            for z in range(dims[2]):
                v = pack(x, y, z)
                if x + 1 < dims[0]:
                    es.append((v, v + STEP[0]))
                if y + 1 < dims[1]:
                    es.append((v, v + STEP[1]))
                if z + 1 < dims[2]:
                    es.append((v, v + STEP[2]))
    return es


# ----------------------------------------------------------------------
# vertex-retirement transfer DP: exact parity-classified polynomials
# ----------------------------------------------------------------------
def parity_polynomial(dims, odd_targets: frozenset[int], K: int) -> list[int]:
    """[v^j] count of edge subsets of the box whose odd-degree vertex set is
    exactly odd_targets, truncated at order K.  Vertex-retirement DP."""
    edges = box_edges(dims)
    last_touch: dict[int, int] = {}
    for idx, (u, w) in enumerate(edges):
        last_touch[u] = idx
        last_touch[w] = idx
    retire_after: dict[int, list[int]] = {}
    for v, idx in last_touch.items():
        retire_after.setdefault(idx, []).append(v)
    # isolated vertices (1x1x1 box) never appear in any edge: they must not
    # be odd targets
    touched = set(last_touch)
    for v in box_vertices(dims):
        if v not in touched and v in odd_targets:
            return [0] * (K + 1)

    states: dict[frozenset[int], list[int]] = {frozenset(): [1] + [0] * K}
    for idx, (u, w) in enumerate(edges):
        new: dict[frozenset[int], list[int]] = {}

        def add(mask, poly):
            cur = new.get(mask)
            if cur is None:
                new[mask] = poly[:]
            else:
                for j, c in enumerate(poly):
                    cur[j] += c

        for mask, poly in states.items():
            add(mask, poly)  # skip the edge
            shifted = [0] + poly[:-1]  # take the edge: one more power of v
            add(mask ^ frozenset((u, w)), shifted)
        for v in retire_after.get(idx, ()):  # retire finished vertices
            want_odd = v in odd_targets
            filtered: dict[frozenset[int], list[int]] = {}
            for mask, poly in new.items():
                if (v in mask) == want_odd:
                    filtered[mask - {v}] = poly
            new = filtered
        states = new
    assert list(states) in ([frozenset()], []), "unretired vertices remain"
    return states.get(frozenset(), [0] * (K + 1))


def poly_div_trunc(num, den, K: int) -> list[int]:
    assert den[0] == 1
    out = [0] * (K + 1)
    for n in range(K + 1):
        s = num[n] if n < len(num) else 0
        for j in range(1, n + 1):
            if j < len(den) and den[j]:
                s -= den[j] * out[n - j]
        out[n] = s
    return out


def box_two_point(dims, x0: int, x1: int, K: int) -> list[int]:
    num = parity_polynomial(dims, frozenset((x0, x1)), K)
    den = parity_polynomial(dims, frozenset(), K)
    return poly_div_trunc(num, den, K)


def box_chi(dims, origin: int, K: int) -> list[int]:
    total = [0] * (K + 1)
    total[0] = 1
    for x in box_vertices(dims):
        if x == origin:
            continue
        tp = box_two_point(dims, origin, x, K)
        for j, c in enumerate(tp):
            total[j] += c
    return total


def brute_two_point(dims, x0: int, x1: int):
    """Complete enumeration of all edge subsets; validation reference."""
    es = box_edges(dims)
    E = len(es)
    num = [0] * (E + 1)
    den = [0] * (E + 1)
    target = frozenset((x0, x1))
    for maskbits in range(1 << E):
        par: dict[int, int] = {}
        m = maskbits
        k = 0
        cnt = 0
        while m:
            if m & 1:
                a, b = es[k]
                cnt += 1
                par[a] = par.get(a, 0) ^ 1
                par[b] = par.get(b, 0) ^ 1
            m >>= 1
            k += 1
        odd = frozenset(v for v, p in par.items() if p)
        if not odd:
            den[cnt] += 1
        elif odd == target:
            num[cnt] += 1
    return num, den


# ----------------------------------------------------------------------
# order-4 dressing obstruction: exact plaquette counts
# ----------------------------------------------------------------------
def plaquettes_meeting(U) -> int:
    seen: set[frozenset[int]] = set()
    axes = ((0, 1), (0, 2), (1, 2))
    for u in U:
        for a, b in axes:
            for da in (0, -1):
                for db in (0, -1):
                    v = u + da * STEP[a] + db * STEP[b]
                    seen.add(
                        frozenset(
                            (v, v + STEP[a], v + STEP[b], v + STEP[a] + STEP[b])
                        )
                    )
    return len(seen)


# ----------------------------------------------------------------------
# inequality audits on an exact coefficient list
# ----------------------------------------------------------------------
def audit(a: list[int]) -> dict[str, object]:
    N = len(a) - 1
    sub_viol = []
    sup_holds = []
    min_margin = None
    for m in range(1, N + 1):
        for n in range(m, N + 1 - m):
            lhs = a[m + n]
            rhs = a[m] * a[n]
            if lhs > rhs:
                sub_viol.append((m, n))
            if lhs >= rhs:
                sup_holds.append((m, n))
            margin = rhs - lhs
            if min_margin is None or margin < min_margin[0]:
                min_margin = (margin, m, n)
    logc_fail = [
        n for n in range(1, N) if a[n] * a[n] < a[n - 1] * a[n + 1]
    ]
    ratios = [Fraction(a[n + 1], a[n]) for n in range(N)]
    ratio_increase = [n for n in range(len(ratios) - 1) if ratios[n + 1] > ratios[n]]
    ratio_over_6 = [n for n in range(len(ratios)) if ratios[n] > 6]
    anchored_super_fail = [
        k for k in range(1, N // 2 + 1) if a[2 * k] < a[k] * a[k]
    ]
    return {
        "order": N,
        "submultiplicativity_violations": sub_viol,
        "supermultiplicativity_pairs_holding": sup_holds,
        "smallest_supermult_counterexample": {
            "pair": [1, 1],
            "a_2": str(a[2]),
            "a_1_squared": str(a[1] * a[1]),
            "margin": str(a[1] * a[1] - a[2]),
        },
        "min_submult_margin": {
            "value": str(min_margin[0]),
            "pair": [min_margin[1], min_margin[2]],
        },
        "log_concavity_failures": logc_fail,
        "ratio_increase_indices": ratio_increase,
        "ratio_over_six_indices": ratio_over_6,
        "anchored_supermult_failures_k": anchored_super_fail,
        "ratio_min": str(min(ratios)),
        "ratio_max": str(max(ratios)),
    }


# ----------------------------------------------------------------------
def main() -> None:
    checks: list[dict[str, object]] = []
    t_start = time.process_time()

    def record(name: str, passed: bool, detail: str) -> None:
        checks.append({"name": name, "passed": bool(passed), "detail": detail})
        print(("PASS" if passed else "FAIL"), name + ":", detail)
        if not passed:
            raise AssertionError(f"check failed: {name}")

    payload = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
    a_own = [int(s) for s in payload["data"]["coefficients"]["a"]]
    a_ext = [int(s) for s in payload["data"]["external"]["oeis_A002913"]]
    N_own = len(a_own) - 1
    record(
        "input_prefix_consistency",
        a_ext[: N_own + 1] == a_own,
        f"theorem-grade a_0..a_{N_own} from e218 agree with the [EXTERNAL] "
        "table prefix",
    )

    # -- audits ------------------------------------------------------------
    own = audit(a_own)
    ext = audit(a_ext)
    record(
        "submultiplicativity_all_pairs",
        not own["submultiplicativity_violations"]
        and not ext["submultiplicativity_violations"],
        f"a_(m+n) <= a_m a_n at every checkable pair: theorem-grade m+n <= "
        f"{N_own}, [EXTERNAL]-grade m+n <= {len(a_ext)-1}; finite evidence "
        "only, promoted to nothing",
    )
    record(
        "supermultiplicativity_smallest_counterexample",
        a_own[2] == 30 and a_own[1] ** 2 == 36 and not own["supermultiplicativity_pairs_holding"],
        "a_2 = 30 < 36 = a_1^2 (theorem-grade, hand-countable: 30 = 2*15 "
        "wedge classes, 6 = 2*3 edge classes); no pair (m,n) with "
        "m+n <= 10 satisfies a_(m+n) >= a_m a_n",
    )
    record(
        "supermultiplicativity_dead_external_range",
        not ext["supermultiplicativity_pairs_holding"],
        "[EXTERNAL] no pair with m+n <= 32 satisfies supermultiplicativity "
        "either",
    )
    record(
        "anchored_supermult_all_anchors_fail",
        ext["anchored_supermult_failures_k"] == list(range(1, 17))
        and own["anchored_supermult_failures_k"] == list(range(1, 6)),
        "a_(2k) < a_k^2 for every testable anchor k (theorem-grade k <= 5, "
        "[EXTERNAL] k <= 16): every anchored Fekete supermultiplicativity "
        "route is refuted at its smallest self-pair",
    )
    record(
        "log_concavity_checked_range",
        not own["log_concavity_failures"] and not ext["log_concavity_failures"],
        "a_n^2 >= a_(n-1) a_(n+1) at every checkable n (equality at n=2); "
        "[CONJECTURE] beyond the checked range",
    )
    record(
        "ratio_monotone_and_ceiling",
        not own["ratio_increase_indices"]
        and not ext["ratio_increase_indices"]
        and not ext["ratio_over_six_indices"],
        "ratios a_(n+1)/a_n are nonincreasing from 6 = a_1/a_0; the unique "
        "checked plateau is a_2/a_1 = a_3/a_2 = 5",
    )

    # -- box study -----------------------------------------------------------
    t0 = time.process_time()
    x0 = pack(0, 0, 0)
    e1 = pack(1, 0, 0)
    # DP validation against complete brute force
    for dims in ((2, 2, 1), (2, 2, 2)):
        for x1 in box_vertices(dims):
            if x1 == x0:
                continue
            num, den = brute_two_point(dims, x0, x1)
            ref = poly_div_trunc(num, den, K_BOX)
            got = box_two_point(dims, x0, x1, K_BOX)
            assert ref == got, (dims, x1)
    record(
        "dp_matches_brute_force",
        True,
        "vertex-retirement DP == complete edge-subset enumeration on "
        "2x2x1 and 2x2x2, every target, through v^8",
    )

    boxes = [(2, 1, 1), (2, 2, 1), (2, 2, 2), (3, 2, 2), (3, 3, 2), (3, 3, 3)]
    tables: dict[str, dict[str, list[str]]] = {}
    for dims in boxes:
        corner = box_chi(dims, x0, K_BOX)
        cx = ((dims[0] - 1) // 2, (dims[1] - 1) // 2, (dims[2] - 1) // 2)
        center = box_chi(dims, pack(*cx), K_BOX)
        tables["x".join(map(str, dims))] = {
            "corner_origin": [str(v) for v in corner],
            "center_origin": [str(v) for v in center],
            "center_origin_at": list(cx),
        }
    chi_221 = [int(s) for s in tables["2x2x1"]["corner_origin"]]
    record(
        "box_negative_coefficient_counterexample",
        chi_221 == [1, 2, 2, 2, 0, -2, -2, -2, 0],
        "2x2x1 corner box: chi_box = 1 + 2v + 2v^2 + 2v^3 - 2v^5 - 2v^6 "
        "- 2v^7 + O(v^9); [v^5] = -2 < 0 refutes box-coefficient "
        "nonnegativity (hand-checkable: 1 + (2v+2v^2+2v^3)/(1+v^4))",
    )
    tp_211 = box_two_point((2, 1, 1), x0, e1, K_BOX)
    tp_221 = box_two_point((2, 2, 1), x0, e1, K_BOX)
    record(
        "box_monotonicity_counterexample",
        tp_211[5] == 0 and tp_221[5] == -1,
        "[v^5]<s_0 s_e1> drops from 0 (2x1x1) to -1 (2x2x1) although "
        "2x1x1 is a subgraph sharing both sites: coefficient-wise "
        "volume monotonicity (the naive GKS bridge premise) is false",
    )
    dom_ok = True
    worst = Fraction(0)
    for name, tab in tables.items():
        for key in ("corner_origin", "center_origin"):
            for n, s in enumerate(tab[key]):
                vabs = abs(int(s))
                if vabs > a_own[n]:
                    dom_ok = False
                if a_own[n]:
                    worst = max(worst, Fraction(vabs, a_own[n]))
    record(
        "box_absolute_domination",
        dom_ok,
        f"|[v^n] chi_box| <= a_n for every tested box, origin, and n <= "
        f"{K_BOX}; worst ratio {worst} = "
        f"{float(worst):.6f} (premise P1' of the conditional bridge lemma: "
        "finite evidence only)",
    )
    stages_box = time.process_time() - t0

    # -- order-4 dressing obstruction ---------------------------------------
    U1 = frozenset((x0, e1))
    U2_adj = frozenset((e1, pack(2, 0, 0)))
    U2_far = frozenset((pack(10, 0, 0), pack(11, 0, 0)))
    p_single = plaquettes_meeting(U1)
    p_adjacent = plaquettes_meeting(U1 | U2_adj)
    p_far = plaquettes_meeting(U1 | U2_far)
    record(
        "dressing_factorization_obstruction",
        p_single == 20 and p_adjacent == 28 and p_far == 40,
        "order-4 dressing: 20 plaquettes meet one edge; a concatenated "
        "adjacent pair is dressed by 28 != 20+20, while a distance-10 pair "
        "is dressed by exactly 40 = 20+20: dressing factorizes at long "
        "distance and fails at contact, which blocks every SAW-style "
        "concatenation injection at coefficient level",
    )

    elapsed = time.process_time() - t_start
    peak_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    record(
        "resource_budget",
        elapsed < 600.0 and peak_rss < 2 * 1024**3,
        f"process time {elapsed:.1f}s < 600s, peak RSS "
        f"{peak_rss/1024**2:.0f} MiB < 2048 MiB",
    )

    payload["meta"]["e219"] = {
        "script": SCRIPT,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "interpreter": ".venv/bin/python",
        "arithmetic": "[COMPUTATION] exact Python integers and Fractions",
        "elapsed_process_seconds": elapsed,
        "peak_rss_bytes": peak_rss,
        "box_stage_seconds": stages_box,
    }
    payload["data"]["inequalities"] = {
        "status": "[COMPUTATION] exhaustive pair audits; [THEOREM] the "
        "smallest supermultiplicativity counterexample (1,1); every "
        "positive direction remains [CONJECTURE] beyond the checked range",
        "theorem_grade": own,
        "external_grade": ext,
    }
    payload["data"]["boxes"] = {
        "status": "[COMPUTATION] exact box susceptibility coefficients "
        "(vertex-retirement DP, brute-force validated)",
        "truncation_order": K_BOX,
        "chi_tables": tables,
        "negative_coefficient_witness": {
            "box": "2x2x1",
            "origin": "corner",
            "order": 5,
            "value": -2,
        },
        "monotonicity_witness": {
            "pair": ["2x1x1", "2x2x1"],
            "observable": "[v^5] <s_0 s_(1,0,0)>",
            "values": [0, -1],
        },
    }
    payload["data"]["obstruction"] = {
        "status": "[COMPUTATION] exact plaquette dressing counts",
        "single_edge": p_single,
        "adjacent_concatenation": p_adjacent,
        "distance_10_pair": p_far,
    }
    fresh = {check["name"] for check in checks}
    payload["checks"] = [
        check for check in payload["checks"] if check["name"] not in fresh
    ] + checks
    RESULT_PATH.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    print(f"extended {RESULT_PATH.relative_to(ROOT)}; {elapsed:.1f} CPU s")


if __name__ == "__main__":
    main()
