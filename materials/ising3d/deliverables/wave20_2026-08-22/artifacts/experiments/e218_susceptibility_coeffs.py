r"""e218 -- Exact simple-cubic Ising susceptibility HT coefficients from finite graphs.

Derives the zero-field susceptibility high-temperature coefficients

    chi(v) = sum_x <s_0 s_x>_{Z^3} = sum_n a_n v^n,   v = tanh K,

exactly, order by order, from finite connected edge clusters, with no input
from any published series.  The identity used is proved in
proofs/susceptibility_coeffs.md:

  [LEMMA cluster identity]  On every finite graph Lambda,
      <s_0 s_x>_Lambda = sum_{C} v^{|C|} * g_Lambda(V(C)),
  where C runs over connected edge sets with odd-degree vertex set exactly
  {0, x}, and g_Lambda(U) = Z^{Lambda \ U} / Z^{Lambda} with
  Z^{X} = sum over fully even edge subsets of E(X) of v^{|edges|}.

  [LEMMA g recursion]  1/g(U) = 1 + sum_{B} v^{|B|} g'(V(B) \ U), where B runs
  over nonempty even edge sets all of whose components meet U, and g' lives on
  the graph with U removed.

  [LEMMA locality/stabilization]  Modulo v^{n+1} every quantity is a finite
  sum over clusters within graph distance n of its base set, so the box
  coefficients stabilize and define the Z^3 coefficients a_n.

Everything is exact Python integer arithmetic.  The order reached is chosen
by a measured process-time calibration (never by a benchmark value), between
N = 9 and N = 10.

Independent internal controls (all recorded as computed checks):
  * full brute-force parity-enumeration identity test on the 2x2x2 box
    (all 7 source-target pairs, through v^8);
  * brute-force identity test on the 3x2x2 box for two targets;
  * hand-countable census values (3 edges, 15 wedges, 20 stars, 3 plaquettes,
    20 four-cycles meeting a fixed edge);
  * an independent depth-first rooted SAW enumeration, with the theorem-grade
    facts a_n = c_n for n <= 4 and a_n < c_n for 5 <= n <= N;
  * the N = 8 census pass is required to be a prefix of the final census.

External comparison (clearly tagged [EXTERNAL], not used in derivation):
OEIS A002913 b-file values (Butera--Comi; Campostrini--Pelissetto--Rossi--
Vicari; Fujiwara--Arisue for n >= 26).

Budgets: single process, peak RSS well under 2 GiB, census-stage budget
1800 CPU s and total producer budget 3000 CPU s.
"""

from __future__ import annotations

import json
import resource
import time
from datetime import datetime, timezone
from pathlib import Path

SCRIPT = "experiments/e218_susceptibility_coeffs.py"
ROOT = Path(__file__).resolve().parents[1]
RESULT_PATH = ROOT / "results" / "bounds" / "susceptibility_coeffs.json"

CENSUS_BUDGET_S = 1800.0
NODE_RATIO_9_OVER_8 = 7.5   # conservative measured growth factors
NODE_RATIO_10_OVER_8 = 60.0

# [EXTERNAL] OEIS A002913 (b-file, retrieved 2026-08-22): HT susceptibility
# coefficients for the spin-1/2 simple-cubic Ising model in v = tanh K.
# n <= 25: P. Butera, M. Comi, Phys. Rev. B 62 (2000) 14837 (through v^23) and
# M. Campostrini, A. Pelissetto, P. Rossi, E. Vicari, Phys. Rev. E 65 (2002)
# 066127 (through v^25); n = 26..32: T. Fujiwara and H. Arisue, workshop
# slides (weaker provenance grade; unused in any theorem-grade claim).
EXTERNAL_A002913 = [
    1, 6, 30, 150, 726, 3510, 16710, 79494, 375174, 1769686,
    8306862, 38975286, 182265822, 852063558, 3973784886, 18527532310,
    86228667894, 401225368086, 1864308847838, 8660961643254,
    40190947325670, 186475398518726, 864404776466406, 4006394107568934,
    18554916271112254, 85923704942057238, 397637244058624494,
    1839992653230056950, 8509528288325589438, 39350934581190850230,
    181885145332015353030, 840628109226856546326, 3883554493872938687622,
]

# [EXTERNAL] OEIS A001412 (already used by proofs/kc_bounds.md): rooted
# simple-cubic SAW counts, comparison values for the internal SAW DFS.
EXTERNAL_A001412 = [
    1, 6, 30, 150, 726, 3534, 16926, 81390, 387966, 1853886, 8809878,
]

# ----------------------------------------------------------------------
# packed geometry: vertex (x, y, z) with |coordinate| <= 12 packs into
# ((x+12)<<10) | ((y+12)<<5) | (z+12); edge = (min_vertex << 2) | axis.
# ----------------------------------------------------------------------
OFF = 12
STEP = (1 << 10, 1 << 5, 1)


def pack(x: int, y: int, z: int) -> int:
    return ((x + OFF) << 10) | ((y + OFF) << 5) | (z + OFF)


def unpack(v: int) -> tuple[int, int, int]:
    return ((v >> 10) - OFF, ((v >> 5) & 31) - OFF, (v & 31) - OFF)


def edge_of(v: int, axis: int) -> int:
    return (v << 2) | axis


def edge_endpoints(e: int) -> tuple[int, int]:
    v = e >> 2
    return v, v + STEP[e & 3]


def edges_at(v: int) -> list[int]:
    out = []
    for axis in range(3):
        out.append((v << 2) | axis)
        out.append(((v - STEP[axis]) << 2) | axis)
    return out


def edge_neighbors(e: int) -> set[int]:
    u, w = edge_endpoints(e)
    out: set[int] = set()
    for vv in (u, w):
        for f in edges_at(vv):
            if f != e:
                out.add(f)
    return out


# ----------------------------------------------------------------------
# census: connected edge clusters up to translation.  The anchor is the
# lexicographically least edge (in packed integer order), fixed at the
# origin edge (0, axis); growth only ever adds edges greater than the
# anchor, so every translation class is generated exactly once.
# ----------------------------------------------------------------------
def enumerate_anchored(nmax: int, visit) -> None:
    for axis in range(3):
        root = edge_of(pack(0, 0, 0), axis)
        u0, w0 = edge_endpoints(root)
        parity = {u0: 1, w0: 1}
        state = {"odd": 2}
        cluster = [root]
        seen = {root}

        def grow(cand: list[int]) -> None:
            visit(cluster, state["odd"], parity)
            if len(cluster) >= nmax:
                return
            for i in range(len(cand)):
                e = cand[i]
                a, b = edge_endpoints(e)
                for vv in (a, b):
                    p = parity.get(vv, 0) ^ 1
                    parity[vv] = p
                    state["odd"] += 1 if p else -1
                cluster.append(e)
                new = [f for f in edge_neighbors(e) if f > root and f not in seen]
                seen.update(new)
                grow(cand[i + 1:] + new)
                seen.difference_update(new)
                cluster.pop()
                for vv in (a, b):
                    p = parity[vv] ^ 1
                    parity[vv] = p
                    state["odd"] += 1 if p else -1

        cand0 = sorted(f for f in edge_neighbors(root) if f > root)
        seen.update(cand0)
        grow(cand0)


# ----------------------------------------------------------------------
# even clusters meeting a vertex set (the correction alphabet)
# ----------------------------------------------------------------------
def plaquettes_meeting(U, forbidden, inside=None) -> int:
    """Exact number of 4-cycles with a vertex in U, no vertex forbidden."""
    seen: set[frozenset[int]] = set()
    axes = ((0, 1), (0, 2), (1, 2))
    for u in U:
        for a, b in axes:
            for da in (0, -1):
                for db in (0, -1):
                    v = u + da * STEP[a] + db * STEP[b]
                    quad = (v, v + STEP[a], v + STEP[b], v + STEP[a] + STEP[b])
                    ok = True
                    for q in quad:
                        if q in forbidden or (inside is not None and not inside(q)):
                            ok = False
                            break
                    if ok:
                        seen.add(frozenset(quad))
    return len(seen)


def even_clusters_meeting(U, budget: int, forbidden, inside=None):
    """Connected even edge clusters D with 4 <= |D| <= budget, V(D) meeting U,
    V(D) disjoint from forbidden, V(D) inside the allowed region.  Each D is
    produced exactly once (anchored at its least U-incident edge)."""
    out = []
    if budget < 4:
        return out

    def vertex_ok(v: int) -> bool:
        if v in forbidden:
            return False
        return inside is None or inside(v)

    def e_ok(f: int) -> bool:
        a, b = edge_endpoints(f)
        return vertex_ok(a) and vertex_ok(b)

    roots = sorted({e for u in U for e in edges_at(u) if e_ok(e)})
    for ridx, root in enumerate(roots):
        banned = set(roots[:ridx])
        u0, w0 = edge_endpoints(root)
        parity = {u0: 1, w0: 1}
        state = {"odd": 2}
        vcount = {u0: 1, w0: 1}
        cluster = [root]
        seen = {root}

        def grow(cand: list[int]) -> None:
            m = len(cluster)
            if state["odd"] == 0 and m >= 4:
                out.append((m, frozenset(vcount)))
            if m >= budget:
                return
            for i in range(len(cand)):
                e = cand[i]
                a, b = edge_endpoints(e)
                for vv in (a, b):
                    p = parity.get(vv, 0) ^ 1
                    parity[vv] = p
                    state["odd"] += 1 if p else -1
                    vcount[vv] = vcount.get(vv, 0) + 1
                cluster.append(e)
                if len(cluster) + state["odd"] // 2 <= budget:
                    new = [
                        f
                        for f in edge_neighbors(e)
                        if f not in seen and f not in banned and e_ok(f)
                    ]
                    seen.update(new)
                    grow(cand[i + 1:] + new)
                    seen.difference_update(new)
                cluster.pop()
                for vv in (a, b):
                    p = parity[vv] ^ 1
                    parity[vv] = p
                    state["odd"] += 1 if p else -1
                    c = vcount[vv] - 1
                    if c:
                        vcount[vv] = c
                    else:
                        del vcount[vv]

        cand0 = [
            f
            for f in sorted(edge_neighbors(root))
            if f not in seen and f not in banned and e_ok(f)
        ]
        seen.update(cand0)
        grow(cand0)
    return out


# ----------------------------------------------------------------------
# g series: g(U) = Z^{X \ U} / Z^X mod v^{K+1}, exact integers
# ----------------------------------------------------------------------
def poly_inv_one_plus(q: list[int], K: int) -> list[int]:
    assert q[0] == 0
    out = [0] * (K + 1)
    out[0] = 1
    for n in range(1, K + 1):
        s = 0
        for j in range(1, n + 1):
            if j < len(q) and q[j]:
                s += q[j] * out[n - j]
        out[n] = -s
    return out


def g_series(U, K: int, forbidden=frozenset(), inside=None, memo=None,
             cross_check=False) -> list[int]:
    if memo is None:
        memo = {}
    key = (U, forbidden, K)
    if key in memo:
        return memo[key]
    if K < 4:
        res = [1] + [0] * K
        memo[key] = res
        return res
    if K < 6:
        # only single 4-cycles can contribute; inner factors are trivial
        p4 = plaquettes_meeting(U, forbidden, inside)
        if cross_check:
            engine = even_clusters_meeting(U, K, forbidden, inside)
            assert len(engine) == p4 and all(t == 4 for t, _ in engine), (
                "plaquette fast path disagrees with the generic engine"
            )
        res = [1] + [0] * K
        res[4] = -p4
        memo[key] = res
        return res
    assert K < 12, "pair-only assembly is valid for K < 12 only"
    comps = even_clusters_meeting(U, K, forbidden, inside)
    q = [0] * (K + 1)
    F2 = frozenset(forbidden | U)
    for idx, (t, vs) in enumerate(comps):
        inner = g_series(frozenset(vs - U), K - t, F2, inside, memo)
        for j, cj in enumerate(inner):
            if t + j <= K and cj:
                q[t + j] += cj
        for jdx in range(idx + 1, len(comps)):
            t2, vs2 = comps[jdx]
            if t + t2 <= K and not (vs & vs2):
                inner2 = g_series(frozenset((vs | vs2) - U), K - t - t2, F2,
                                  inside, memo)
                for j, cj in enumerate(inner2):
                    if t + t2 + j <= K and cj:
                        q[t + t2 + j] += cj
    res = poly_inv_one_plus(q, K)
    memo[key] = res
    return res


# ----------------------------------------------------------------------
# finite-box validation machinery
# ----------------------------------------------------------------------
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
                    es.append(edge_of(v, 0))
                if y + 1 < dims[1]:
                    es.append(edge_of(v, 1))
                if z + 1 < dims[2]:
                    es.append(edge_of(v, 2))
    return es


def brute_two_point(dims, x0: int, x1: int):
    """Exact numerator/denominator of <s_x0 s_x1> on the free-boundary box,
    by complete enumeration of edge subsets classified by odd-vertex set."""
    es = box_edges(dims)
    E = len(es)
    endpoints = [edge_endpoints(e) for e in es]
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
                a, b = endpoints[k]
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


def cluster_two_point(dims, x0: int, x1: int, K: int) -> list[int]:
    """Cluster-identity evaluation of <s_x0 s_x1> on the box, mod v^{K+1}."""
    vset = set(box_vertices(dims))
    inside = vset.__contains__
    es = set(box_edges(dims))
    all_edges = sorted(es)
    total = [0] * (K + 1)
    memo: dict = {}
    target = tuple(sorted((x0, x1)))
    for root in all_edges:
        u0, w0 = edge_endpoints(root)
        parity = {u0: 1, w0: 1}
        state = {"odd": 2}
        cluster = [root]
        seen = {root}

        def visit() -> None:
            if state["odd"] != 2:
                return
            odds = tuple(sorted(v for v, p in parity.items() if p))
            if odds != target:
                return
            m = len(cluster)
            verts = frozenset(v for e in cluster for v in edge_endpoints(e))
            gs = g_series(verts, K - m, frozenset(), inside, memo)
            for j, c in enumerate(gs):
                if m + j <= K:
                    total[m + j] += c

        def grow(cand: list[int]) -> None:
            visit()
            if len(cluster) >= K:
                return
            for i in range(len(cand)):
                e = cand[i]
                a, b = edge_endpoints(e)
                for vv in (a, b):
                    p = parity.get(vv, 0) ^ 1
                    parity[vv] = p
                    state["odd"] += 1 if p else -1
                cluster.append(e)
                new = [
                    f
                    for f in edge_neighbors(e)
                    if f in es and f > root and f not in seen
                ]
                seen.update(new)
                grow(cand[i + 1:] + new)
                seen.difference_update(new)
                cluster.pop()
                for vv in (a, b):
                    p = parity[vv] ^ 1
                    parity[vv] = p
                    state["odd"] += 1 if p else -1

        cand0 = [f for f in sorted(edge_neighbors(root)) if f in es and f > root]
        seen.update(cand0)
        grow(cand0)
    return total


# ----------------------------------------------------------------------
# infinite-lattice census + assembly
# ----------------------------------------------------------------------
def run_census(N: int):
    """Return (histogram[(size, odd)] -> classes, two_odd_by_size,
    small = [(size, vertexset)] for size <= N-4)."""
    hist: dict[tuple[int, int], int] = {}
    two_odd = [0] * (N + 1)
    small: list[tuple[int, frozenset[int]]] = []

    def vis(cluster, odd, parity):
        m = len(cluster)
        key = (m, odd)
        hist[key] = hist.get(key, 0) + 1
        if odd == 2:
            two_odd[m] += 1
            if m <= N - 4:
                small.append(
                    (m, frozenset(v for e in cluster for v in edge_endpoints(e)))
                )

    enumerate_anchored(N, vis)
    return hist, two_odd, small


def assemble_coefficients(N: int, two_odd, small, cross_check_fastpath: bool):
    a = [0] * (N + 1)
    a[0] = 1
    for m in range(1, N + 1):
        a[m] += 2 * two_odd[m]
    memo: dict = {}
    for m, verts in small:
        gs = g_series(verts, N - m, frozenset(), None, memo,
                      cross_check=cross_check_fastpath and (N - m) <= 5)
        for j in range(4, N - m + 1):
            if gs[j]:
                a[m + j] += 2 * gs[j]
    return a


def saw_counts(N: int) -> list[int]:
    """Rooted SAW counts c_0..c_N on Z^3 by direct visited-set DFS."""
    counts = [0] * (N + 1)
    counts[0] = 1
    origin = pack(0, 0, 0)
    visited = {origin}
    deltas = [STEP[0], STEP[1], STEP[2], -STEP[0], -STEP[1], -STEP[2]]

    def rec(v: int, depth: int) -> None:
        if depth == N:
            return
        for d in deltas:
            w = v + d
            if w in visited:
                continue
            counts[depth + 1] += 1
            visited.add(w)
            rec(w, depth + 1)
            visited.remove(w)

    rec(origin, 0)
    return counts


# ----------------------------------------------------------------------
def main() -> None:
    checks: list[dict[str, object]] = []
    stages: dict[str, float] = {}

    def record(name: str, passed: bool, detail: str) -> None:
        checks.append({"name": name, "passed": bool(passed), "detail": detail})
        print(("PASS" if passed else "FAIL"), name + ":", detail)
        if not passed:
            raise AssertionError(f"check failed: {name}")

    t_start = time.process_time()

    # -- V1: hand-countable alphabet facts ------------------------------
    t0 = time.process_time()
    U_edge = frozenset(edge_endpoints(edge_of(pack(0, 0, 0), 0)))
    p4_edge = plaquettes_meeting(U_edge, frozenset())
    record(
        "four_cycles_meeting_edge",
        p4_edge == 20,
        f"4-cycles meeting a fixed edge: {p4_edge} == 12+12-4",
    )
    engine4 = even_clusters_meeting(U_edge, 4, frozenset())
    record(
        "engine_matches_plaquette_count",
        len(engine4) == 20 and all(t == 4 for t, _ in engine4),
        f"generic even-cluster engine at budget 4 returns {len(engine4)} 4-cycles",
    )
    stages["alphabet_controls"] = time.process_time() - t0

    # -- V2: full brute-force identity on 2x2x2 -------------------------
    t0 = time.process_time()
    dims = (2, 2, 2)
    x0 = pack(0, 0, 0)
    K_val = 8
    all_ok = True
    for x1 in box_vertices(dims):
        if x1 == x0:
            continue
        num, den = brute_two_point(dims, x0, x1)
        ref = poly_div_trunc(num, den, K_val)
        got = cluster_two_point(dims, x0, x1, K_val)
        if ref != got:
            all_ok = False
            break
    record(
        "identity_2x2x2_all_targets",
        all_ok,
        "cluster identity == brute-force parity enumeration on 2x2x2, "
        "all 7 targets, through v^8",
    )
    stages["identity_2x2x2"] = time.process_time() - t0

    # -- V3: brute-force identity on 3x2x2, two targets -----------------
    t0 = time.process_time()
    dims = (3, 2, 2)
    ok = True
    for x1 in (pack(1, 0, 0), pack(2, 1, 1)):
        num, den = brute_two_point(dims, x0, x1)
        ref = poly_div_trunc(num, den, K_val)
        got = cluster_two_point(dims, x0, x1, K_val)
        if ref != got:
            ok = False
            break
    record(
        "identity_3x2x2_two_targets",
        ok,
        "cluster identity == brute-force parity enumeration on 3x2x2 "
        "for x=(1,0,0) and x=(2,1,1), through v^8",
    )
    stages["identity_3x2x2"] = time.process_time() - t0

    # -- calibration census at N = 8 ------------------------------------
    t0 = time.process_time()
    hist8, two_odd8, _ = run_census(8)
    t_census8 = time.process_time() - t0
    stages["census_8"] = t_census8
    record(
        "census_hand_counts",
        hist8[(1, 2)] == 3
        and hist8[(2, 2)] == 15
        and hist8[(3, 4)] == 20
        and hist8[(3, 2)] == 75
        and hist8[(4, 0)] == 3,
        "3 edge classes, 15 wedges = C(6,2), 20 stars = C(6,3), "
        "75 two-odd 3-clusters, 3 plaquette classes",
    )

    projected_10 = t_census8 * NODE_RATIO_10_OVER_8
    projected_9 = t_census8 * NODE_RATIO_9_OVER_8
    if projected_10 <= CENSUS_BUDGET_S:
        N = 10
    elif projected_9 <= CENSUS_BUDGET_S:
        N = 9
    else:
        N = 8
    print(
        f"calibration: census(8) = {t_census8:.2f}s, projections "
        f"9 -> {projected_9:.1f}s, 10 -> {projected_10:.1f}s, chosen N = {N}"
    )

    # -- main census -----------------------------------------------------
    t0 = time.process_time()
    if N == 8:
        hist, two_odd, small = hist8, two_odd8, None
        raise SystemExit("census budget too small for a useful run")
    hist, two_odd, small = run_census(N)
    stages["census_main"] = time.process_time() - t0
    record(
        "census_prefix_regression",
        all(hist.get(k, 0) == v for k, v in hist8.items() if k[0] <= 8)
        and two_odd[:9] == two_odd8[:9],
        f"N={N} census agrees with the independent N=8 pass on every "
        "size <= 8 class count",
    )

    # -- corrections and assembly ----------------------------------------
    t0 = time.process_time()
    a = assemble_coefficients(N, two_odd, small, cross_check_fastpath=True)
    stages["corrections"] = time.process_time() - t0

    record(
        "leading_orders_hand",
        a[0] == 1 and a[1] == 6 and a[2] == 30,
        "a_0 = 1; a_1 = 2*3 = 6 (three directed edge classes, two roots); "
        "a_2 = 2*15 = 30 (fifteen wedge classes)",
    )
    even_g_only = all(
        all(g == 0 for j, g in enumerate(gs) if j % 2 == 1)
        for gs in [g_series(frozenset(edge_endpoints(edge_of(pack(0, 0, 0), d))), N - 1)
                   for d in range(3)]
    )
    record(
        "bipartite_parity_structure",
        even_g_only,
        "g-series of an edge has zero odd-order coefficients "
        "(bipartite lattice: even clusters have even size)",
    )

    # -- SAW comparison ----------------------------------------------------
    t0 = time.process_time()
    c = saw_counts(N)
    stages["saw_dfs"] = time.process_time() - t0
    record(
        "saw_dfs_vs_A001412",
        c == EXTERNAL_A001412[: N + 1],
        f"[EXTERNAL cross-check] internal SAW DFS c_0..c_{N} equals OEIS "
        "A001412 prefix",
    )
    record(
        "a_equals_c_through_4",
        all(a[n] == c[n] for n in range(5)),
        "a_n = c_n for n <= 4 (first even cluster needs 4 edges, so "
        "corrections start at order 5)",
    )
    record(
        "a_below_c_from_5",
        all(a[n] < c[n] for n in range(5, N + 1)),
        f"a_n < c_n for 5 <= n <= {N}: exact margins "
        f"{[c[n] - a[n] for n in range(5, N + 1)]}",
    )
    record(
        "a5_correction_identity",
        a[5] == 2 * two_odd[5] - 2 * 3 * 20
        and 2 * two_odd[5] == c[5] + 2 * 48,
        "a_5 = 2*1815 - 120 = 3510: 1815 = 1767 path classes + 48 "
        "plaquette-plus-pendant classes; correction = -2*(3 edge classes)*"
        "(20 four-cycles each)",
    )

    # -- external comparison ----------------------------------------------
    record(
        "external_A002913_prefix",
        a == EXTERNAL_A002913[: N + 1],
        f"[EXTERNAL cross-check] derived a_0..a_{N} equals the published "
        "series prefix (OEIS A002913 b-file; Butera--Comi / Campostrini "
        "et al. values); the derivation used no published input",
    )

    elapsed = time.process_time() - t_start
    peak_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    record(
        "resource_budget",
        elapsed < 3000.0 and peak_rss < 2 * 1024**3,
        f"process time {elapsed:.1f}s < 3000s, peak RSS "
        f"{peak_rss/1024**2:.0f} MiB < 2048 MiB",
    )

    payload = {
        "meta": {
            "script": SCRIPT,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": ".venv/bin/python",
            "arithmetic": "[COMPUTATION] exact Python integers throughout; "
            "no floating point enters any stored value",
            "elapsed_process_seconds": elapsed,
            "peak_rss_bytes": peak_rss,
            "stages_process": stages,
            "order_reached": N,
            "census_budget_seconds": CENSUS_BUDGET_S,
        },
        "data": {
            "coefficients": {
                "status": "[COMPUTATION] exact, derived from finite graphs only",
                "convention": "chi(v) = sum_x <s_0 s_x> on Z^3, v = tanh K, "
                "free-boundary stabilized coefficients "
                "(proofs/susceptibility_coeffs.md, Lemma 1-3)",
                "a": [str(x) for x in a],
                "order": N,
            },
            "census": {
                "status": "[COMPUTATION]",
                "classes_by_size_and_oddcount": {
                    f"{m},{o}": hist[(m, o)] for (m, o) in sorted(hist)
                },
                "two_odd_by_size": two_odd,
            },
            "saw": {
                "status": "[COMPUTATION] internal DFS; "
                "[EXTERNAL] OEIS A001412 comparison",
                "c": [str(x) for x in c],
                "first_strict_gap": 5,
                "gaps_c_minus_a": [str(c[n] - a[n]) for n in range(N + 1)],
            },
            "external": {
                "status": "[EXTERNAL] comparison only, unused in derivation",
                "oeis_A002913": [str(x) for x in EXTERNAL_A002913],
                "provenance": "n <= 25 peer-reviewed (Butera--Comi PRB 62 "
                "(2000) 14837; Campostrini--Pelissetto--Rossi--Vicari PRE 65 "
                "(2002) 066127); n = 26..32 Fujiwara--Arisue workshop slides "
                "(weaker grade, never used in theorem-grade claims)",
            },
        },
        "checks": checks,
    }
    RESULT_PATH.parent.mkdir(parents=True, exist_ok=True)
    RESULT_PATH.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    print(f"wrote {RESULT_PATH.relative_to(ROOT)}; N = {N}, "
          f"{elapsed:.1f} CPU s, peak {peak_rss/1024**2:.0f} MiB")


if __name__ == "__main__":
    main()
