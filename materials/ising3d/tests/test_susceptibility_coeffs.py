"""Independent exact verifier for results/bounds/susceptibility_coeffs.json.

Imports no producer module.  All geometry uses plain coordinate tuples (the
producers use packed integers), and every re-derivation uses a different
algorithm than the producer that created the corresponding stored value:

  * census of connected edge-cluster translation classes, sizes <= 6, by a
    canonical-form breadth-first closure over classes (producer: anchored
    untried-set depth-first enumeration);
  * census sizes 7 and 8 by an independently written recursive anchored
    enumeration over tuple coordinates;
  * susceptibility coefficients a_0..a_8 re-assembled from the census with
    closed-form order-4/order-6 corrections (only 4-cycle and 6-cycle counts
    enter through v^8; the producer used a general recursion);
  * finite-box susceptibility tables by literal spin-sum enumeration
    (producers: parity-subset enumeration and a vertex-retirement DP);
  * inequality audits, completions, envelope, and directed enclosures by
    fresh Fraction arithmetic (bisection integer roots, a shorter atanh
    series with exact tail).

Verification scope: a_0..a_8 are re-derived here in full independence;
a_9 and a_10 are covered by the producer's internal identity controls plus
exact agreement with the [EXTERNAL] published table, and are re-checked here
only against that table.  Run time is a few minutes.
"""

from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "bounds" / "susceptibility_coeffs.json"

FAILURES: list[str] = []

EXTERNAL_A002913 = [
    1, 6, 30, 150, 726, 3510, 16710, 79494, 375174, 1769686,
    8306862, 38975286, 182265822, 852063558, 3973784886, 18527532310,
    86228667894, 401225368086, 1864308847838, 8660961643254,
    40190947325670, 186475398518726, 864404776466406, 4006394107568934,
    18554916271112254, 85923704942057238, 397637244058624494,
    1839992653230056950, 8509528288325589438, 39350934581190850230,
    181885145332015353030, 840628109226856546326, 3883554493872938687622,
]
INCUMBENT_LOWER = Fraction(2122159753270231627267517174278577806020, 10**40)


def check(name: str, passed: bool, detail: str = "") -> None:
    print(f"{'PASS' if passed else 'FAIL'} {name} -- {detail}")
    if not passed:
        FAILURES.append(name)


# ----------------------------------------------------------------------
# tuple geometry
# ----------------------------------------------------------------------
AXES = ((1, 0, 0), (0, 1, 0), (0, 0, 1))


def eplus(v, d):
    return (v[0] + d[0], v[1] + d[1], v[2] + d[2])


def eminus(v, d):
    return (v[0] - d[0], v[1] - d[1], v[2] - d[2])


def edge_endpoints(e):
    v, a = e
    return v, eplus(v, AXES[a])


def edges_at(v):
    out = []
    for a in range(3):
        out.append((v, a))
        out.append((eminus(v, AXES[a]), a))
    return out


def edge_neighbors(e):
    u, w = edge_endpoints(e)
    res = set()
    for vv in (u, w):
        for f in edges_at(vv):
            if f != e:
                res.add(f)
    return res


def cluster_vertices(cluster):
    out = set()
    for e in cluster:
        u, w = edge_endpoints(e)
        out.add(u)
        out.add(w)
    return out


def odd_vertices(cluster):
    par = {}
    for e in cluster:
        for vv in edge_endpoints(e):
            par[vv] = par.get(vv, 0) ^ 1
    return {v for v, p in par.items() if p}


def canonical(cluster):
    verts = cluster_vertices(cluster)
    xm = min(v[0] for v in verts)
    ym = min(v[1] for v in verts)
    zm = min(v[2] for v in verts)
    return frozenset(((v[0] - xm, v[1] - ym, v[2] - zm), a) for (v, a) in cluster)


# ----------------------------------------------------------------------
# census algorithm 1: canonical-class breadth-first closure, sizes <= 6
# ----------------------------------------------------------------------
def class_bfs_census(nmax):
    classes = [set() for _ in range(nmax + 1)]
    seeds = {canonical(frozenset([(((0, 0, 0)), a)])) for a in range(3)}
    classes[1] = set(seeds)
    for m in range(1, nmax):
        for cl in classes[m]:
            grown = set()
            for e in cl:
                for f in edge_neighbors(e):
                    if f not in cl:
                        grown.add(f)
            for f in grown:
                classes[m + 1].add(canonical(cl | {f}))
    return classes


# ----------------------------------------------------------------------
# census algorithm 2: recursive anchored enumeration, sizes <= 8
# ----------------------------------------------------------------------
def anchored_census(nmax):
    """histogram[(size, odd_count)] -> class count; also returns the two-odd
    vertex sets for sizes <= 4 (needed by the closed-form corrections)."""
    hist: dict[tuple[int, int], int] = {}
    small: list[tuple[int, frozenset]] = []

    def edge_key(e):
        return (e[0], e[1])  # tuples compare lexicographically

    for axis in range(3):
        root = ((0, 0, 0), axis)
        rk = edge_key(root)
        parity: dict[tuple, int] = {}
        for vv in edge_endpoints(root):
            parity[vv] = 1
        cluster = [root]
        seen = {root}

        def visit():
            odd = sum(parity.values())
            m = len(cluster)
            hist[(m, odd)] = hist.get((m, odd), 0) + 1
            if odd == 2 and m <= 4:
                small.append((m, frozenset(cluster_vertices(cluster))))

        def grow(cand):
            visit()
            if len(cluster) >= nmax:
                return
            for i in range(len(cand)):
                e = cand[i]
                for vv in edge_endpoints(e):
                    parity[vv] = parity.get(vv, 0) ^ 1
                cluster.append(e)
                new = [
                    f
                    for f in edge_neighbors(e)
                    if edge_key(f) > rk and f not in seen
                ]
                seen.update(new)
                grow(cand[i + 1:] + new)
                seen.difference_update(new)
                cluster.pop()
                for vv in edge_endpoints(e):
                    parity[vv] = parity.get(vv, 0) ^ 1

        cand0 = sorted(
            (f for f in edge_neighbors(root) if edge_key(f) > rk), key=edge_key
        )
        seen.update(cand0)
        grow(cand0)
    return hist, small


# ----------------------------------------------------------------------
# closed-form correction alphabet: 4-cycles and 6-cycles meeting a set
# ----------------------------------------------------------------------
def four_cycles_meeting(U):
    seen = set()
    pairs = ((0, 1), (0, 2), (1, 2))
    for u in U:
        for a, b in pairs:
            da, db = AXES[a], AXES[b]
            for sa in (0, -1):
                for sb in (0, -1):
                    v = (u[0] + sa * da[0] + sb * db[0],
                         u[1] + sa * da[1] + sb * db[1],
                         u[2] + sa * da[2] + sb * db[2])
                    quad = (v, eplus(v, da), eplus(v, db), eplus(eplus(v, da), db))
                    seen.add(frozenset(quad))
    return len(seen)


def six_cycles_meeting(U):
    """6-cycles (as edge sets) with a vertex in U, by self-avoiding closed
    walks of length six from each vertex of U."""
    dirs = [AXES[0], AXES[1], AXES[2],
            tuple(-c for c in AXES[0]), tuple(-c for c in AXES[1]),
            tuple(-c for c in AXES[2])]

    def edge_between(u, w):
        for a in range(3):
            if eplus(u, AXES[a]) == w:
                return (u, a)
            if eplus(w, AXES[a]) == u:
                return (w, a)
        raise AssertionError("not adjacent")

    found = set()
    for u in U:
        path = [u]

        def walk():
            if len(path) == 7:
                return
            for d in dirs:
                w = eplus(path[-1], d)
                if len(path) == 6:
                    if w == u:
                        cyc = frozenset(
                            edge_between(path[i], path[(i + 1) % 6])
                            for i in range(6)
                        )
                        if len(cyc) == 6:
                            found.add(cyc)
                    continue
                if w in path or w == u and len(path) < 6:
                    continue
                path.append(w)
                walk()
                path.pop()

        walk()
    return len(found)


# ----------------------------------------------------------------------
# finite boxes by literal spin sums
# ----------------------------------------------------------------------
def box_sites(dims):
    return [
        (x, y, z)
        for x in range(dims[0])
        for y in range(dims[1])
        for z in range(dims[2])
    ]


def box_bonds(dims):
    bonds = []
    sites = box_sites(dims)
    index = {s: i for i, s in enumerate(sites)}
    for s in sites:
        for a in range(3):
            t = eplus(s, AXES[a])
            if t in index:
                bonds.append((index[s], index[t]))
    return sites, index, bonds


def _spin_polys(dims):
    """Per-configuration expansion of prod_bonds (1 + v s_i s_j): returns
    (n_sites, index, list over configurations of (spins, poly))."""
    sites, index, bonds = box_bonds(dims)
    n = len(sites)
    rows = []
    for conf in range(1 << n):
        spins = [1 if (conf >> i) & 1 else -1 for i in range(n)]
        poly = [1]
        for (i, j) in bonds:
            sij = spins[i] * spins[j]
            newp = [0] * (len(poly) + 1)
            for k, c in enumerate(poly):
                newp[k] += c
                newp[k + 1] += c * sij
            poly = newp
        rows.append((spins, poly))
    return n, index, rows


def _ratio_trunc(num, den, K):
    """num/den mod v^{K+1}; den[0] divides everything exactly."""
    scale = den[0]
    out = []
    for k in range(K + 1):
        s = Fraction(num[k] if k < len(num) else 0, scale)
        for j in range(1, k + 1):
            if j < len(den):
                s -= Fraction(den[j], scale) * out[k - j]
        out.append(s)
    assert all(f.denominator == 1 for f in out)
    return [int(f) for f in out]


def spin_two_point(dims, s0, s1, K):
    """<sigma_s0 sigma_s1> mod v^{K+1} by literal enumeration of spins."""
    n, index, rows = _spin_polys(dims)
    i0, i1 = index[s0], index[s1]
    E = len(rows[0][1]) - 1
    num = [0] * (E + 1)
    den = [0] * (E + 1)
    for spins, poly in rows:
        w = spins[i0] * spins[i1]
        for k, c in enumerate(poly):
            den[k] += c
            num[k] += c * w
    assert den[0] == 1 << n
    return _ratio_trunc(num, den, K)


def spin_chi(dims, origin, K):
    """[v^k] sum_x <s_origin s_x> mod v^{K+1}, one spin enumeration."""
    n, index, rows = _spin_polys(dims)
    i0 = index[origin]
    E = len(rows[0][1]) - 1
    den = [0] * (E + 1)
    nums = {i: [0] * (E + 1) for i in range(n) if i != i0}
    for spins, poly in rows:
        s0 = spins[i0]
        for k, c in enumerate(poly):
            den[k] += c
        for i in nums:
            w = s0 * spins[i]
            row = nums[i]
            for k, c in enumerate(poly):
                row[k] += c * w
    assert den[0] == 1 << n
    total = [0] * (K + 1)
    total[0] = 1
    for i, num in nums.items():
        tp = _ratio_trunc(num, den, K)
        for j, c in enumerate(tp):
            total[j] += c
    return total


# ----------------------------------------------------------------------
# rational analysis helpers (independent of the producers')
# ----------------------------------------------------------------------
def introot_bisect(x, n):
    lo, hi = 0, 1
    while hi**n <= x:
        hi <<= 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if mid**n <= x:
            lo = mid
        else:
            hi = mid
    return lo


def atanh_iv(x: Fraction, terms: int = 90):
    total = Fraction(0)
    xsq = x * x
    p = x
    for k in range(terms + 1):
        total += p / (2 * k + 1)
        p *= xsq
    tail = p / ((2 * terms + 3) * (1 - xsq))
    return total, total + tail


def dec_down(f: Fraction, places: int) -> str:
    s = 10**places
    q = (f.numerator * s) // f.denominator
    return f"{q // s}.{q % s:0{places}d}"


def dec_up(f: Fraction, places: int) -> str:
    s = 10**places
    q = -((-f.numerator * s) // f.denominator)
    return f"{q // s}.{q % s:0{places}d}"


# ----------------------------------------------------------------------
def main() -> int:
    payload = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    data = payload["data"]
    a = [int(s) for s in data["coefficients"]["a"]]
    N = data["coefficients"]["order"]
    check("artifact_shape",
          set(payload) >= {"meta", "data", "checks"}
          and len(a) == N + 1 and N >= 9,
          f"meta/data/checks present, order {N}")
    check("all_stored_checks_passed",
          all(c["passed"] for c in payload["checks"]),
          f"{len(payload['checks'])} stored checks")
    check("external_table_matches",
          [int(s) for s in data["external"]["oeis_A002913"]] == EXTERNAL_A002913
          and a == EXTERNAL_A002913[: N + 1],
          "stored external table equals the verifier's copy; stored "
          "coefficients equal its prefix")

    # -- census, two independent algorithms -----------------------------
    classes = class_bfs_census(6)
    bfs_hist = {}
    for m in range(1, 7):
        for cl in classes[m]:
            o = len(odd_vertices(cl))
            bfs_hist[(m, o)] = bfs_hist.get((m, o), 0) + 1
    hist, small4 = anchored_census(8)
    check("census_algorithms_agree",
          all(hist.get(k, 0) == v for k, v in bfs_hist.items())
          and all(bfs_hist.get(k, 0) == v for k, v in hist.items()
                  if k[0] <= 6),
          "class-BFS census == anchored census on every size <= 6 cell")
    hand = {(1, 2): 3, (2, 2): 15, (3, 2): 75, (3, 4): 20, (4, 0): 3,
            (4, 2): 363, (4, 4): 315}
    check("census_hand_values",
          all(hist.get(k, 0) == v for k, v in hand.items()),
          "3 edges, 15 wedges, 20 stars, 3 plaquettes, 75/363 two-odd")
    stored_census = data["census"]["classes_by_size_and_oddcount"]
    check("stored_census_prefix",
          all(stored_census.get(f"{m},{o}", 0) == hist.get((m, o), 0)
              for m in range(1, 9) for o in range(0, 10, 2)),
          "stored census equals the independent census on every size <= 8 "
          "cell")

    # -- coefficients a_0..a_8 from scratch ------------------------------
    two_odd = {m: hist.get((m, 2), 0) for m in range(1, 9)}
    b = [0] * 9
    b[0] = 1
    for m in range(1, 9):
        b[m] = 2 * two_odd[m]
    # corrections through v^8 need only 4-cycle and 6-cycle counts:
    #   order m+4: -2 * #(4-cycles meeting V(C)) for each two-odd C, |C| = m
    #   order m+6: -2 * #(6-cycles meeting V(C))
    for m, verts in small4:
        p4 = four_cycles_meeting(verts)
        if m + 4 <= 8:
            b[m + 4] -= 2 * p4
        if m + 6 <= 8:
            b[m + 6] -= 2 * six_cycles_meeting(verts)
    check("coefficients_rederived_0_to_8",
          b == a[:9],
          f"independent re-derivation gives {b}")
    check("hand_identity_a6",
          two_odd[6] == 8775
          and all(four_cycles_meeting(v) == 28 for m, v in small4 if m == 2)
          and b[6] == 2 * 8775 - 2 * 15 * 28 == 16710,
          "a_6 = 2*8775 - 2*15*28: every wedge is dressed by exactly 28 "
          "four-cycles")
    edge_set = frozenset(((0, 0, 0), (1, 0, 0)))
    check("hand_identity_a5",
          four_cycles_meeting(edge_set) == 20
          and b[5] == 2 * two_odd[5] - 2 * 3 * 20 == 3510,
          "a_5 = 2*1815 - 120")
    check("six_cycles_edge_count",
          six_cycles_meeting(edge_set) == 214,
          "214 six-cycles meet a fixed edge (enters a_7)")

    # -- boxes by literal spin sums ---------------------------------------
    K = data["boxes"]["truncation_order"]
    tables = data["boxes"]["chi_tables"]
    ok = True
    for dims in ((2, 1, 1), (2, 2, 1), (2, 2, 2), (3, 2, 2)):
        name = "x".join(map(str, dims))
        stored = [int(s) for s in tables[name]["corner_origin"]]
        got = spin_chi(dims, (0, 0, 0), K)
        if stored != got:
            ok = False
        cx = tuple((d - 1) // 2 for d in dims)
        stored_c = [int(s) for s in tables[name]["center_origin"]]
        got_c = spin_chi(dims, cx, K)
        if stored_c != got_c:
            ok = False
    check("box_tables_spin_sum",
          ok,
          "literal spin-sum enumeration reproduces the stored chi tables "
          "on 2x1x1, 2x2x1, 2x2x2, 3x2x2 (both origins, through v^8)")
    check("negative_coefficient_witness",
          spin_chi((2, 2, 1), (0, 0, 0), 8) == [1, 2, 2, 2, 0, -2, -2, -2, 0],
          "[v^5] chi(2x2x1, corner) = -2 < 0")
    check("monotonicity_witness",
          spin_two_point((2, 1, 1), (0, 0, 0), (1, 0, 0), 8)[5] == 0
          and spin_two_point((2, 2, 1), (0, 0, 0), (1, 0, 0), 8)[5] == -1,
          "[v^5]<s_0 s_e1> drops 0 -> -1 from 2x1x1 to 2x2x1")
    dom = True
    for name, tab in tables.items():
        for key in ("corner_origin", "center_origin"):
            for n, s in enumerate(tab[key]):
                if abs(int(s)) > a[n]:
                    dom = False
    check("absolute_domination_restated",
          dom, "|[v^n] chi_box| <= a_n on every stored table")
    obstruction = data["obstruction"]
    adj = frozenset(((0, 0, 0), (1, 0, 0), (2, 0, 0)))
    far = frozenset(((0, 0, 0), (1, 0, 0), (10, 0, 0), (11, 0, 0)))
    check("dressing_obstruction_reproduced",
          four_cycles_meeting(edge_set) == obstruction["single_edge"] == 20
          and four_cycles_meeting(adj) == obstruction["adjacent_concatenation"] == 28
          and four_cycles_meeting(far) == obstruction["distance_10_pair"] == 40,
          "order-4 dressing: 20 / 28 (adjacent, != 40) / 40 (distance 10)")
    stored_two_odd = data["census"]["two_odd_by_size"]
    check("stored_two_odd_prefix",
          all(stored_two_odd[m] == two_odd[m] for m in range(1, 9))
          and stored_two_odd[0] == 0,
          "stored two-odd class counts equal the independent census "
          "for every size <= 8")

    # -- inequality audits -------------------------------------------------
    ext = EXTERNAL_A002913
    Ne = len(ext) - 1
    sub_viol = [(m, n) for m in range(1, Ne + 1) for n in range(m, Ne + 1 - m)
                if ext[m + n] > ext[m] * ext[n]]
    sup_hold = [(m, n) for m in range(1, Ne + 1) for n in range(m, Ne + 1 - m)
                if ext[m + n] >= ext[m] * ext[n]]
    check("submultiplicativity_reaudit",
          not sub_viol, "no violation with m+n <= 32")
    check("supermultiplicativity_reaudit",
          not sup_hold and ext[2] == 30 and ext[1] ** 2 == 36,
          "fails at every pair; smallest counterexample (1,1): 30 < 36")
    check("log_concavity_reaudit",
          all(ext[n] * ext[n] >= ext[n - 1] * ext[n + 1]
              for n in range(1, Ne)),
          "log-concave on the checked range")
    ratios = [Fraction(ext[n + 1], ext[n]) for n in range(Ne)]
    check("ratio_reaudit",
          all(ratios[n + 1] <= ratios[n] for n in range(Ne - 1))
          and [n for n in range(Ne - 1) if ratios[n + 1] == ratios[n]] == [1]
          and max(ratios) == 6,
          "ratios nonincreasing from 6; unique plateau 5,5")
    check("anchored_supermult_reaudit",
          all(ext[2 * k] < ext[k] ** 2 for k in range(1, 17)),
          "a_(2k) < a_k^2 for k <= 16")
    ineq = data["inequalities"]
    check("stored_inequality_tables",
          not ineq["external_grade"]["submultiplicativity_violations"]
          and not ineq["external_grade"]["supermultiplicativity_pairs_holding"]
          and ineq["external_grade"]["smallest_supermult_counterexample"]["margin"] == "6"
          and not ineq["theorem_grade"]["submultiplicativity_violations"],
          "stored audit tables match the recomputation")

    # -- directed section ---------------------------------------------------
    directed = data["directed"]
    check("no_bound_claimed",
          directed["claimed_new_bound"] is None,
          "no all-n inequality is proved, so no new K_c bound is claimed")
    targets = directed["conditional_route"]["targets"]
    ok = True
    for ns, tab in targets.items():
        n = int(ns)
        q = (10**40) ** n // ext[n]
        r = introot_bisect(q, n)
        blo, bhi = Fraction(r, 10**40), Fraction(r + 1, 10**40)
        if [str(blo), str(bhi)] != tab["bracket_v"]:
            ok = False
        alo, _ = atanh_iv(blo)
        _, ahi = atanh_iv(bhi)
        lo40, hi40 = tab["atanh_enclosure_40dp"]
        if not (Fraction(lo40) <= alo and ahi <= Fraction(hi40)
                and dec_down(alo, 40) == lo40 and dec_up(ahi, 40) == hi40):
            ok = False
    check("directed_targets_reproduced",
          ok, f"brackets and 40-decimal enclosures re-derived for n in "
              f"{sorted(int(k) for k in targets)}")
    beats = {}
    for n in range(1, Ne + 1):
        q = (10**40) ** n // ext[n]
        r = introot_bisect(q, n)
        alo, _ = atanh_iv(Fraction(r, 10**40))
        beats[n] = alo > INCUMBENT_LOWER
    check("crossing_order_reproduced",
          min(n for n in beats if beats[n]) == 21
          and directed["conditional_route"]["first_order_beating_incumbent"] == 21
          and not beats[N]
          and directed["conditional_route"]["beats_incumbent_by_order"]
          == {str(n): bool(beats[n]) for n in range(1, Ne + 1)},
          "conditional route beats the incumbent from n = 21 only, "
          "never within the theorem-grade range")
    alo6, ahi6 = atanh_iv(Fraction(1, 6))
    check("ratio_ceiling_valueless_reproduced",
          ahi6 < INCUMBENT_LOWER,
          f"atanh(1/6) <= {dec_up(ahi6, 20)} < incumbent")

    # -- completions ---------------------------------------------------------
    H = directed["completions"]["closure_horizon"]
    s = list(ext)
    for n in range(Ne + 1, H + 1):
        s.append(min(s[m] * s[n - m] for m in range(1, n)))
    check("closure_reproduced",
          [str(x) for x in s[33:41]] == directed["completions"]["closure_tail_33_to_40"]
          and all(s[m + n] <= s[m] * s[n]
                  for m in range(1, H) for n in range(m, H + 1 - m)),
          f"min-plus closure re-derived and fully submultiplicative to {H}")
    d = list(ext) + [6**n for n in range(Ne + 1, H + 1)]
    first_viol = next((m, n) for m in range(1, H)
                      for n in range(m, H + 1 - m) if d[m + n] > d[m] * d[n])
    check("divergent_completion_reproduced",
          list(first_viol) == directed["completions"]["divergent_first_violation_pair"]
          and d[: Ne + 1] == ext,
          f"first violation at {first_viol}, invisible below order 33")
    check("envelope_reproduced",
          all(ext[n] <= (1000 * (n + 1) ** 3) ** (n + 2)
              for n in range(1, Ne + 1))
          and all(6**n <= (1000 * (n + 1) ** 3) ** (n + 2)
                  for n in range(1, H + 1))
          and all(s[n] <= (1000 * (n + 1) ** 3) ** (n + 2)
                  for n in range(1, H + 1))
          and data["directed"]["completions"]["envelope"]
          == "E_n = (1000 (n+1)^3)^(n+2) (Lemma 6)",
          "all values and both completions respect the provable envelope "
          "E_n = (1000 (n+1)^3)^(n+2)")

    print()
    print("scope note: a_0..a_8 re-derived fully independently here; "
          "a_9, a_10 rest on the producer's brute-force identity controls "
          "plus exact agreement with the published table.")
    if FAILURES:
        print(f"FAIL ({len(FAILURES)}): {', '.join(FAILURES)}")
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
