"""Standalone independent audit of the e140 upper-beyond deliverables.

Everything decisive is recomputed by a route that differs from the experiment:

* the polycube census through n = 10 cells (complete for contour areas <= 28
  by Loomis-Whitney) is regenerated with a different canonical encoding, a
  different growth routine, and its own per-class complement flood fill; the
  wave-6 table (areas <= 24) is re-derived and compared verbatim with
  results/bounds/peierls_upper.json, and the tree strata are re-derived both
  from the bond strata and by an independent leaf-extension enumeration
  (through n = 10);
* every exponential enclosure e^{-2K} is recomputed by the RECIPROCAL route
  1/e^{2K} with an all-positive Taylor series and an exact geometric remainder
  bound (the producer used the alternating-series route); the I_3 anchor is
  recomputed through mpmath interval arithmetic and nested against the
  repo-certified interval of results/bounds/upper_infrared.json;
* the limitation inequality, the K* and K_26 bisections, the lambda-tail rows,
  and the certified 1e-6-grid head+tail endpoint are re-decided on exact
  Fraction endpoints of the test's own enclosures;
* every row of the artifact's universal tables (level counts, bond strata,
  exact N(A), tree classes, tail components including the n = 10 rows, checks)
  is asserted individually against the test's own regenerated census data.

Run: PYTHONPATH=src .venv/bin/python tests/test_upper_beyond.py
"""

from __future__ import annotations

import json
import resource
import sys
import time
from fractions import Fraction
from pathlib import Path

# the artifact carries exact Fractions whose numerators reach ~10^4 digits
sys.set_int_max_str_digits(1_000_000)

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "bounds" / "upper_beyond.json"
WAVE6 = ROOT / "results" / "bounds" / "peierls_upper.json"
INFRARED = ROOT / "results" / "bounds" / "upper_infrared.json"

N_MAX = 10          # full census depth in this test (complete for areas <= 28)
TREE_N_MAX = 10     # independent tree census depth
TAIL_AREAS = (28, 30, 32, 34, 36, 38, 40, 42)
HEAD_AREAS = (6, 10, 14, 16, 18, 20, 22, 24, 26)
BUDGET_CENSUS_S = 3600.0
BUDGET_TREE_S = 2400.0
EXP_ORDER = 50

# fixed polycubes (translations quotiented, orientations distinct): OEIS A001931
EXPECTED_LEVELS = {1: 1, 2: 3, 3: 15, 4: 86, 5: 534, 6: 3481, 7: 23502, 8: 162913, 9: 1152870, 10: 8294738}
WAVE6_HEAD = {6: 1, 8: 0, 10: 6, 12: 0, 14: 45, 16: 12, 18: 332, 20: 240, 22: 2538, 24: 3040}
EXPECTED_TREES = {1: 1, 2: 3, 3: 15, 4: 83, 5: 486, 6: 2967, 7: 18748, 8: 121725, 9: 807381, 10: 5447203}
N26 = 20448
N28 = 33240
WAVE6_K0 = Fraction(848449, 500000)
LAMBDA = Fraction(11**11, 10**10)
DIRS = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))

CAVITY_WITNESS = [
    (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1),
    (1, 1, 0), (-1, 1, 0), (-1, -1, 0), (0, -1, 1), (0, -1, -1),
]


def cpu() -> float:
    return time.process_time()


def rss_gib() -> float:
    raw = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return (raw if sys.platform == "darwin" else raw * 1024) / 1024**3


# ---------------------------------------------------------------------------
# Independent census engine: 24-bit packed coordinates, full re-canonicalisation
# ---------------------------------------------------------------------------

def canon(cells) -> bytes:
    """Translate minima to zero, sort, pack each cell as three bytes."""
    xs = [c[0] for c in cells]
    ys = [c[1] for c in cells]
    zs = [c[2] for c in cells]
    mx, my, mz = min(xs), min(ys), min(zs)
    keys = sorted(((c[0] - mx) << 16) | ((c[1] - my) << 8) | (c[2] - mz) for c in cells)
    out = bytearray()
    for k in keys:
        out.append((k >> 16) & 255)
        out.append((k >> 8) & 255)
        out.append(k & 255)
    return bytes(out)


def unpack(data: bytes):
    return [(data[3 * i], data[3 * i + 1], data[3 * i + 2]) for i in range(len(data) // 3)]


def grow(prev: set[bytes]) -> set[bytes]:
    """All canonical (n+1)-cell supersets, by full re-canonicalisation."""
    out: set[bytes] = set()
    add = out.add
    for data in prev:
        cells = unpack(data)
        for (x, y, z) in cells:
            for dx, dy, dz in DIRS:
                cand = cells + [(x + dx, y + dy, z + dz)]
                if len(set(cand)) == len(cells) + 1:
                    add(canon(cand))
    return out


def bond_count(cells) -> int:
    """Number of shared faces, counting the three positive directions."""
    s = set(cells)
    return sum(
        1
        for (x, y, z) in s
        for d in ((1, 0, 0), (0, 1, 0), (0, 0, 1))
        if (x + d[0], y + d[1], z + d[2]) in s
    )


def complement_connected(cells) -> bool:
    """Flood fill of the complement inside the one-cell-expanded box.

    Flat-array BFS with index order (z-major, then y, then x) -- a different
    memory layout from the experiment's (x-major) flood fill.
    """
    cells = set(cells)
    xs = [c[0] for c in cells]
    ys = [c[1] for c in cells]
    zs = [c[2] for c in cells]
    lo = (min(xs) - 1, min(ys) - 1, min(zs) - 1)
    Dx = max(xs) - lo[0] + 2
    Dy = max(ys) - lo[1] + 2
    Dz = max(zs) - lo[2] + 2
    grid = bytearray(Dz * Dy * Dx)
    for (x, y, z) in cells:
        grid[((z - lo[2]) * Dy + (y - lo[1])) * Dx + (x - lo[0])] = 1
    XY = Dy * Dx
    mark = bytearray(len(grid))
    mark[0] = 1
    stack = [0]
    reached = 1
    while stack:
        p = stack.pop()
        z, r = divmod(p, XY)
        y, x = divmod(r, Dx)
        for q in (
            (p - XY if z > 0 else -1), (p + XY if z < Dz - 1 else -1),
            (p - Dx if y > 0 else -1), (p + Dx if y < Dy - 1 else -1),
            (p - 1 if x > 0 else -1), (p + 1 if x < Dx - 1 else -1),
        ):
            if q >= 0 and not grid[q] and not mark[q]:
                mark[q] = 1
                stack.append(q)
                reached += 1
    return reached + len(cells) == len(grid)


def census(n_max: int, budget_s: float):
    """Levels 1..n_max with bond strata (all and cavity-free)."""
    t0 = cpu()
    level = {canon([(0, 0, 0)])}
    counts = {1: 1}
    strata_all = {(1, 0): 1}
    strata_ok = {(1, 0): 1}
    n = 1
    while n < n_max:
        n += 1
        level = grow(level)
        counts[n] = len(level)
        for data in level:
            cells = unpack(data)
            b = bond_count(cells)
            strata_all[(n, b)] = strata_all.get((n, b), 0) + 1
            if complement_connected(cells):
                strata_ok[(n, b)] = strata_ok.get((n, b), 0) + 1
        print(f"  test census n={n}: {counts[n]} (cpu {cpu() - t0:.1f}s)", flush=True)
        if cpu() - t0 > budget_s:
            raise SystemExit(f"test census exceeded budget {budget_s}s at n={n}")
    return counts, strata_all, strata_ok


def tree_levels(n_max: int, budget_s: float) -> dict[int, int]:
    """Independent leaf-extension tree census (adjacent to exactly one cell)."""
    t0 = cpu()
    level = {canon([(0, 0, 0)])}
    trees = {1: 1}
    for n in range(1, n_max):
        nxt: set[bytes] = set()
        for data in level:
            cells = unpack(data)
            members = set(cells)
            for (x, y, z) in cells:
                for dx, dy, dz in DIRS:
                    c = (x + dx, y + dy, z + dz)
                    if c in members:
                        continue
                    nbrs = [
                        (c[0] + d[0], c[1] + d[1], c[2] + d[2]) in members for d in DIRS
                    ]
                    if sum(nbrs) == 1:
                        nxt.add(canon(cells + [c]))
        level = nxt
        trees[n + 1] = len(level)
        print(f"  test trees n={n + 1}: {len(level)} (cpu {cpu() - t0:.1f}s)", flush=True)
        if cpu() - t0 > budget_s:
            raise SystemExit(f"test tree census exceeded budget {budget_s}s at n={n + 1}")
    return trees


def binomial(n: int, k: int) -> int:
    if k < 0 or k > n:
        return 0
    r = 1
    for i in range(k):
        r = r * (n - i) // (i + 1)
    return r

# ---------------------------------------------------------------------------
# Independent exponential enclosures: reciprocal of the positive series
# ---------------------------------------------------------------------------

def exp_pos_bounds(u: Fraction, order: int = EXP_ORDER) -> tuple[Fraction, Fraction]:
    """Exact S <= e^{u} <= S + rem for 0 <= u <= 1/2 (all-positive series)."""
    term = Fraction(1)
    s = Fraction(1)
    for j in range(1, order + 1):
        term = term * u / j
        s += term
    # remainder: ratio of successive omitted terms <= u/(order+2) <= 1/104
    first = term * u / (order + 1)
    ratio = u / (order + 2)
    if ratio >= 1:
        raise AssertionError("remainder ratio out of range")
    rem = first / (1 - ratio)
    return s, s + rem


def exp_neg_bounds(t: Fraction, q: int = 10**48) -> tuple[Fraction, Fraction]:
    """Exact lo <= e^{-t} <= hi for rational t >= 0 via 1/e^{t}.

    The tight enclosure is then rounded OUTWARD to the grid 1/q (lo down, hi
    up): the value stays certified inside, and downstream series evaluations
    run on bounded-denominator Fractions instead of gcd-heavy giant ones.
    """
    if t < 0:
        raise ValueError("t >= 0 required")
    if t == 0:
        return Fraction(1), Fraction(1)
    n_sub = max(1, -(-2 * t.numerator // t.denominator))
    u = t / n_sub
    if u > Fraction(1, 2):
        n_sub += 1
        u = t / n_sub
    lo_e, hi_e = exp_pos_bounds(u)
    lo, hi = 1 / hi_e**n_sub, 1 / lo_e**n_sub
    lo_q, hi_q = lo * q, hi * q
    return (
        Fraction(lo_q.numerator // lo_q.denominator, q),
        Fraction(-((-hi_q.numerator) // hi_q.denominator), q),
    )


def binomial(n: int, k: int) -> int:
    if k < 0 or k > n:
        return 0
    r = 1
    for i in range(k):
        r = r * (n - i) // (i + 1)
    return r


def frac(text: str) -> Fraction:
    num, den = text.split("/")
    return Fraction(int(num), int(dec))


def main() -> None:
    t0 = cpu()
    art = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    data = art["data"]

    # ---- producer self-checks all passed, no floats in certificate ----------
    assert all(c["passed"] for c in art["checks"]), [
        c["name"] for c in art["checks"] if not c["passed"]
    ]
    assert art["meta"]["arithmetic"]["floating_point_used_in_certificate"] is False
    assert "ruled out" in art["meta"]["constraint_class"]
    print("artifact checks all passed, exact-arithmetic flags present: PASS")

    # ---- independent census through n = 9 ------------------------------------
    counts, strata_all, strata_ok = census(N_MAX, BUDGET_CENSUS_S)
    for n, expected in EXPECTED_LEVELS.items():
        assert counts[n] == expected, (n, counts[n], expected)
    art_levels = {int(k): v for k, v in data["census"]["level_counts"].items()}
    for n in range(1, N_MAX + 1):
        assert art_levels[n] == counts[n], (n, art_levels[n], counts[n])
    for n in range(1, N_MAX + 1):
        total = sum(v for (nn, _), v in strata_all.items() if nn == n)
        assert total == counts[n]
    assert all(strata_all.get((n, b), 0) == 0 for n in range(1, N_MAX + 1) for b in range(n - 1))
    print("independent census levels 1..10 + strata sanity: PASS")

    # ---- full bond-strata table comparison (every n <= 10 row) ---------------
    art_strata = data["census"]["bond_strata"]
    for key_text, row in art_strata.items():
        n, b = (int(t) for t in key_text.split(","))
        assert row["all"] == strata_all.get((n, b), 0), (key_text, row["all"], strata_all.get((n, b), 0))
        assert row["cavity_free"] == strata_ok.get((n, b), 0), (key_text, row["cavity_free"], strata_ok.get((n, b), 0))
    assert len(art_strata) == len(strata_all), (len(art_strata), len(strata_all))
    print(f"full bond strata table ({len(art_strata)} rows) regenerated: PASS")

    # ---- flood-fill micro-tests and cavity-free strata ------------------------
    assert complement_connected([(x, y, z) for x in (0, 1) for y in (0, 1) for z in (0, 1)])
    plate = [(x, y, 0) for x in range(3) for y in range(3) if (x, y) != (1, 1)]
    assert complement_connected(plate)
    assert not complement_connected(CAVITY_WITNESS)
    assert bond_count(CAVITY_WITNESS) == 10
    n_cav = sum(strata_all[k] - strata_ok.get(k, 0) for k in strata_all)
    assert n_cav == 0, n_cav
    assert data["census"]["all_strata_cavity_free"] is True
    assert data["census"]["min_cavity_search"]["min_connected_cavity_cell_count"] == 11
    print("flood-fill witnesses + cavity-free strata n<=10 + min cavity 11: PASS")

    # ---- exact N(A) table rows, wave-6 verbatim, decompositions --------------
    Ntab: dict[int, int] = {}
    for (n, b), c in strata_ok.items():
        A = 6 * n - 2 * b
        Ntab[A] = Ntab.get(A, 0) + n * c
    for A, v in WAVE6_HEAD.items():
        assert Ntab.get(A, 0) == v, (A, Ntab.get(A, 0), v)
    w6 = json.loads(WAVE6.read_text(encoding="utf-8"))["data"]["exact_head"]["N_of_A"]
    for k, v in w6.items():
        assert Ntab.get(int(k), 0) == v, (k, Ntab.get(int(k), 0), v)
    art_N = {int(k): v for k, v in data["exact_N_table"].items()}
    assert set(art_N) == set(range(6, 29, 2)) - {8, 12}, sorted(art_N)
    for A in range(6, 27, 2):
        assert art_N.get(A, 0) == Ntab.get(A, 0), (A, art_N.get(A), Ntab.get(A))
    assert art_N[28] == N28 == Ntab[28]
    assert Ntab[26] == N26 == 6 * strata_ok[(6, 5)] + 7 * strata_ok[(7, 8)]
    assert N28 == 7 * strata_ok[(7, 7)] + 8 * strata_ok[(8, 10)] + 9 * strata_ok[(9, 13)]
    dec26 = data["N26_decomposition"]
    assert dec26["6*T_6"] == 6 * strata_ok[(6, 5)] and dec26["7*c(7,8)"] == 7 * strata_ok[(7, 8)]
    dec28 = data["N28_decomposition"]
    assert dec28["7*c(7,7)"] == 7 * strata_ok[(7, 7)]
    assert dec28["8*c(8,10)"] == 8 * strata_ok[(8, 10)]
    assert dec28["9*c(9,13)"] == 9 * strata_ok[(9, 13)]
    assert dec28["10*c(10,16)"] == 0
    # Loomis-Whitney completeness rows
    assert 216 * 9**2 <= 26**3 < 216 * 10**2
    assert 216 * 10**2 <= 28**3 < 216 * 11**2
    print("N(A) rows A<=26, wave-6 verbatim, N26/N28 decompositions: PASS")

    # ---- tree classes: identity rows + independent census through n = 10 ------
    trees = {n: strata_ok[(n, n - 1)] for n in range(1, N_MAX + 1)}
    for n in range(1, N_MAX + 1):
        assert 6 * n - 2 * (n - 1) == 4 * n + 2
    art_T = {int(k): v for k, v in data["tree_classes"]["T_n"].items()}
    for n, v in trees.items():
        assert art_T[n] == v == EXPECTED_TREES[n], (n, art_T[n], v)
    trees10 = tree_levels(TREE_N_MAX, BUDGET_TREE_S)
    for n in range(1, TREE_N_MAX + 1):
        assert trees10[n] == art_T[n] == EXPECTED_TREES[n], (n, trees10[n], art_T[n])
    ratios = data["tree_classes"]["successive_ratios"]
    for n in range(1, 10):
        f = Fraction(ratios[str(n)])
        assert f == Fraction(art_T[n + 1], art_T[n])
    print("tree classes n<=10 by independent leaf extension + ratios: PASS")

    # ---- tail components: EVERY row (n <= 10) regenerated by the test -------
    M_art = {}
    for A in TAIL_AREAS:
        row = data["tail_lower_bound_components"][str(A)]
        M_art[A] = row["M_A_n_le_10"]
        total = 0
        for n_text, c in row["strata"].items():
            n = int(n_text)
            b = 3 * n - A // 2
            assert 6 * n - 2 * b == A, (A, n, b)
            # every row, including the n = 10 rows, is regenerated exactly
            assert strata_ok.get((n, b), 0) == c, (A, n, b, c)
            total += n * c
        assert total == M_art[A], (A, total, M_art[A])
    assert M_art[28] == N28 == Ntab[28]
    assert M_art[42] == 10 * art_T[10] == 10 * EXPECTED_TREES[10] == 10 * trees10[10]
    m30_regen = 7 * strata_ok[(7, 6)] + 8 * strata_ok[(8, 9)] + 9 * strata_ok[(9, 12)] + 10 * strata_ok[(10, 15)]
    assert m30_regen == M_art[30]
    print("tail component rows (every row n<=10 regenerated exactly): PASS")

    # ---- exponential enclosures vs mpmath and the repo-certified anchor -------
    import mpmath as mp
    from decimal import Decimal

    mp.iv.dps = 100
    q = mp.iv
    I3_iv = (
        q.sqrt(6)
        * q.gamma(q.mpf(1) / 24)
        * q.gamma(q.mpf(5) / 24)
        * q.gamma(q.mpf(7) / 24)
        * q.gamma(q.mpf(11) / 24)
        / (96 * q.pi**3)
    )
    with mp.workdps(200):

        def ivf(v) -> Fraction:
            f = mp.mpf(v)
            sign, man, exp, _bc = f._mpf_
            out = Fraction(man) * Fraction(2) ** exp
            return -out if sign else out

        i3_lo = ivf(I3_iv.a)
        i3_hi = ivf(I3_iv.b)
    inf = json.loads(INFRARED.read_text(encoding="utf-8"))
    repo_lo = Fraction(Decimal(inf["data"]["certified_constants"]["I3"][0]))
    repo_hi = Fraction(Decimal(inf["data"]["certified_constants"]["I3"][1]))
    assert repo_lo - Fraction(1, 10**90) <= i3_lo and i3_hi <= repo_hi + Fraction(1, 10**90)
    kbar_lo, kbar_hi = i3_lo / 2, i3_hi / 2
    art_kbar = data["limitation_certificate"]["kbar_interval"]
    assert Fraction(art_kbar["lower_exact"]) < kbar_hi and kbar_lo < Fraction(art_kbar["upper_exact"])
    for tt in (Fraction(1), Fraction(1, 2), Fraction(7, 2), 2 * kbar_hi):
        lo, hi = exp_neg_bounds(tt)
        with mp.workdps(80):
            ref = mp.exp(-mp.mpf(tt.numerator) / tt.denominator)
            assert mp.mpf(lo.numerator) / lo.denominator < ref < mp.mpf(hi.numerator) / hi.denominator
        assert hi - lo < Fraction(1, 10**40)
    print("I_3 anchor + reciprocal-route exp enclosures: PASS")

    # ---- x* enclosure and the master limitation inequality --------------------
    x_hi_t = exp_neg_bounds(i3_lo)[1]   # e^{-i3_lo} bounds x* from above
    x_lo_t = exp_neg_bounds(i3_hi)[0]   # e^{-i3_hi} bounds x* from below
    art_x = data["limitation_certificate"]["x_star_interval"]
    ax_lo, ax_hi = Fraction(art_x["lower_exact"]), Fraction(art_x["upper_exact"])
    assert x_lo_t < ax_hi and ax_lo < x_hi_t  # enclosures overlap
    head = {A: Ntab[A] for A in HEAD_AREAS}

    def series(coeffs: dict[int, int], y: Fraction) -> Fraction:
        return sum(Fraction(c) * y**A for A, c in coeffs.items())

    h26_lo = series(head, x_lo_t)
    tail_lo = series(M_art, x_lo_t)
    margin = h26_lo + tail_lo - Fraction(1, 2)
    assert margin > 0, margin
    art_margin_txt = data["limitation_certificate"]["margin_over_one_half_decimal"]
    art_margin_lo = Fraction(art_margin_txt.strip("[]").split(",")[0].strip())
    assert art_margin_lo > 0, art_margin_txt
    assert data["limitation_certificate"]["margin_positive_certified_on_exact_fractions"] is True
    slack_lo = Fraction(1, 2) - h26_lo
    assert tail_lo > slack_lo
    print(f"master inequality reproduced: margin {float(margin):.6f} > 0; tail > slack: PASS")

    # ---- K* bisection with the test's own enclosures --------------------------
    head28 = dict(head)
    head28[28] = N28

    def V_lower(K: Fraction) -> Fraction:
        y, _ = exp_neg_bounds(2 * K)
        m_wo28 = {A: M_art[A] for A in M_art if A != 28}
        return series(head28, y) + series(m_wo28, y)

    assert V_lower(kbar_hi) > Fraction(1, 2)
    lo, hi = kbar_hi, Fraction(3, 2)
    for _ in range(64):
        mid = (lo + hi) / 2
        if V_lower(mid) > Fraction(1, 2):
            lo = mid
        else:
            hi = mid
    art_k = data["limitation_certificate"]["K_star_interval_exact"]
    ak_lo, ak_hi = Fraction(art_k[0]), Fraction(art_k[1])
    assert lo < ak_hi and ak_lo < hi  # brackets overlap
    assert lo - kbar_hi > 0
    assert Fraction(data["limitation_certificate"]["K_star_minus_kbar_exact"]) > 0
    print(f"K* reproduced: {float(lo):.9f} > K_bar = {float(kbar_hi):.9f}: PASS")

    # ---- spurious head-only crossing K_26 and its refutation ------------------
    def H26(K: Fraction, upper: bool) -> Fraction:
        enc = exp_neg_bounds(2 * K)
        y = enc[1] if upper else enc[0]
        return series(head, y)

    b_lo, b_hi = Fraction(1, 8), kbar_lo
    assert H26(b_lo, True) > Fraction(1, 2) > H26(b_hi, False)
    for _ in range(64):
        mid = (b_lo + b_hi) / 2
        if H26(mid, True) < Fraction(1, 2):
            b_hi = mid
        elif H26(mid, False) > Fraction(1, 2):
            b_lo = mid
        else:
            break
    art_b = data["limitation_certificate"]["K26_spurious_interval_exact"]
    assert b_hi > Fraction(art_b[0]) - Fraction(1, 10**9)
    assert b_lo < Fraction(art_b[1]) + Fraction(1, 10**9)
    assert b_hi < kbar_lo  # the feigned gain is real arithmetic, and spurious
    assert V_lower(b_hi) > Fraction(1, 2)
    print(f"K_26 spurious crossing reproduced ({float(b_lo):.9f}) and refuted: PASS")

    # ---- lambda-tail rows and the certified grid endpoint ---------------------
    # Lagrange identity by the test's own series reversion
    deg = 33
    T = [0] * (deg + 1)
    T[1] = 1
    for _ in range(deg + 2):
        base = [1] + T[1:]
        acc = [1] + [0] * deg
        for _p in range(11):
            nxt = [0] * (deg + 1)
            for i, a in enumerate(acc):
                if a:
                    for j, b in enumerate(base):
                        if b and i + j <= deg:
                            nxt[i + j] += a * b
            acc = nxt
        shifted = [0] * (deg + 1)
        for i in range(deg):
            shifted[i + 1] = acc[i]
        if shifted == T:
            break
        T = shifted
    R = [0] * (deg + 1)
    base = [1] + T[1:]
    acc = [1] + [0] * deg
    for _ in range(12):
        nxt = [0] * (deg + 1)
        for i, a in enumerate(acc):
            if a:
                for j, b in enumerate(base):
                    if b and i + j <= deg:
                        nxt[i + j] += a * b
        acc = nxt
    for i in range(deg):
        R[i + 1] = acc[i]
    assert R[1] == 1
    for A in range(2, deg + 1):
        assert R[A] == 12 * binomial(11 * A, A - 2) // (A - 1), A
    for A in range(26, 121, 2):
        assert 100 * binomial(11 * A, A - 2) <= LAMBDA**A, A
        B_A = 12 * binomial(11 * A, A - 2) // (A - 1)
        assert 625 * A * B_A <= 78 * LAMBDA**A, A
    print("Lagrange + binomial + animal rows (own reversion): PASS")

    def cert(K: Fraction, upper: bool) -> Fraction:
        enc = exp_neg_bounds(2 * K)
        y = enc[1] if upper else enc[0]
        ly = LAMBDA * y
        if ly >= 1:
            return Fraction(10**6)
        return series(head28, y) + Fraction(78, 625) * ly**30 / (1 - ly * ly)

    g = Fraction(data["proved_exponential_tail"]["grid_endpoint_exact"])
    assert cert(g, True) < Fraction(1, 2), "grid endpoint must pass (upper)"
    assert cert(g - Fraction(1, 10**6), False) > Fraction(1, 2), "grid-1 must fail (lower)"
    assert g < WAVE6_K0
    conv = data["proved_exponential_tail"]["convergence_threshold_half_log_lambda"]
    conv_lo = Fraction(conv[0])
    assert LAMBDA * exp_neg_bounds(2 * conv_lo)[0] > 1  # tail divergent at conv_lo
    assert conv_lo > kbar_hi
    print(f"grid endpoint {float(g):.6f} re-certified, below wave-6 1.696898: PASS")

    print(f"\nALL CHECKS PASSED (cpu {cpu() - t0:.1f}s, peak RSS {rss_gib():.2f} GiB)")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"FAIL: {exc}")
        raise
