"""Exact 3D outer-contour census through area 28, the certified failure of the
plain Peierls certificate at the infrared endpoint I_3/2, and the exact
threshold structure of the route.

Result classification
---------------------
[COMPUTATION] Exhaustive exact census of rooted outer contours of the simple
cubic Ising model, by cell count n <= 10 (complete for every contour area
A <= 28 by the discrete Loomis-Whitney cutoff), extending the wave-6 exact head
(A <= 24, experiments/e41_peierls_upper.py) by two further areas:

    A          6  8 10 12  14  16  18  20   22   24    26    28
    N(A)       1  0  6  0  45  12 332 240 2538 3040 20448 33240

with the bond strata c(n, b) (translation classes with n cells and b bonds)
fully tabulated for n <= 10.  Every n <= 10 class has a cavity-free complement,
certified DIRECTLY by an exact flood fill of the complement of every one of the
9 638 143 level classes (zero cavity-bearing classes found), and sharpened
by an independent frontier enumeration: the analytic decomposition
|D|=1 (six-spoke front), |D|=2 (disconnected domino neighbourhood), and
3<=|D|<=5 (trio neighbourhood bound 16-|D| >= 11) combined with the direct
computation for |D|>=6 shows the smallest cavity-bearing polycube has exactly
11 cells (a tree witness is exhibited).  The tree strata T_n = c(n, n-1) are
cross-checked by an independent leaf-extension enumeration.

[THEOREM: limitation] Let K_bar = I_3/2 be the certified incumbent upper
endpoint (proofs/upper_infrared.md), x* = e^{-2 K_bar} = e^{-I_3}.  Using ONLY
exact areas A <= 26 plus the n <= 10 components of areas 28, 30, ..., 42
(strict lower bounds M_A on N(A)),

    V(x*) := sum_{A<=26} N(A) x*^A + sum_{A in {28..42}, even} M_A x*^A  >  1/2,

proved as an exact rational inequality between Fraction endpoints of a
100-digit enclosure of x*.  Since every term is nonnegative, the full contour
sum satisfies sum_A N(A) e^{-2KA} >= V(e^{-2K}) > 1/2 for every K <= K_bar, so
the plain Peierls certificate "sum_A N(A) e^{-2KA} < 1/2 implies K_c <= K"
(kc_upper_peierls.md eq. (2.5)-(2.6)) fails at the incumbent and at every
smaller K.  The bisection threshold K* of V satisfies K* > K_bar (exact
rational gap reported): the plain-certificate route loses by a certified
positive margin and cannot reach, let alone beat, I_3/2.  The head-only
crossing K_26 (where the truncated head alone reaches 1/2, deep below K_bar)
is computed exactly and refuted: V(e^{-2 K_26}) > 1/2 by a certified margin.

[COMPUTATION: proved exponential truncation tail where it converges] The
wave-6 universal-cover tail bound N(A) <= (78/625) lambda^A, lambda =
11^11/10^10, is re-derived and machine-verified here (Lagrange-inversion
identity checked by iterative series reversion, binomial majorant checked row
by row), and combined with the extended exact head through A = 28 to certify
the least 10^-6-grid endpoint of this head-plus-tail strategy: an improvement
of the wave-6 number 1.696898, still ~6.7 times the incumbent.  The tail
converges only for K > (1/2) log lambda = 1.6754..., a certified 1.4227 above
K_bar, which quantifies exactly why no finite exact head can help this
particular tail argument (wave-6 Section 5 limitation, restated with the new
head).

The two-point infrared constraint class is ruled out for improvement by
proofs/upper_infrared.md (exact optimum I_3/2, method-optimal) and
proofs/mag_floor.md (zero-magnetization profile saturation, O(1/L) torus LP
floor); this front therefore works strictly in the Peierls contour-count
class and offers no two-point-only variant.

Arithmetic: exact Python integers throughout the census; the only
transcendental inputs are e^{-2K} at rational K (and I_3), enclosed by pure
exact-Fraction Taylor arithmetic (alternating series after subdivision --
independent of any floating-point library), cross-checked against mpmath
interval arithmetic at dps=100 following proofs/upper_infrared.md Section 2.
Every proof-step comparison is decided on exact Fraction endpoints.  No
floating-point value enters any proof step.
"""

from __future__ import annotations

import json
import resource
import sys
import time
from fractions import Fraction
from pathlib import Path

# Exact certificate Fractions on y^A, A <= 42, have numerators up to ~10^6
# digits; only ceilings/divisions, never str(), decide a proof step, but the
# small fraction_text payloads need the stdlib limit lifted for safety.
sys.set_int_max_str_digits(5_000_000)

SCRIPT = "experiments/e140_upper_beyond.py"
ROOT = Path(__file__).resolve().parents[1]
RESULT_PATH = ROOT / "results" / "bounds" / "upper_beyond.json"
WAVE6_PATH = ROOT / "results" / "bounds" / "peierls_upper.json"
INFRARED_PATH = ROOT / "results" / "bounds" / "upper_infrared.json"

N_MAX = 10  # cell cutoff: complete for contour areas A <= 28 (Loomis-Whitney)
TAIL_AREAS = (28, 30, 32, 34, 36, 38, 40, 42)
HEAD_AREAS = (6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26)
IV_DPS = 100
EXP_ORDER = 44  # per-subdivision alternating-Taylor order (u <= 1/2)

# Predeclared stage budgets in CPU seconds (50%-of-CPU machine policy).
BUDGET_CENSUS_S = 3600.0
BUDGET_TREE_S = 2400.0
BUDGET_CAVITY_S = 900.0
BUDGET_ARITH_S = 900.0
RSS_LIMIT_GIB = 5.5

EXPECTED_LEVELS = {
    # fixed polycubes (translations quotiented, orientations distinct): OEIS A001931
    1: 1, 2: 3, 3: 15, 4: 86, 5: 534, 6: 3481, 7: 23502, 8: 162913, 9: 1152870,
    10: 8294738,
}
WAVE6_HEAD = {6: 1, 8: 0, 10: 6, 12: 0, 14: 45, 16: 12, 18: 332, 20: 240, 22: 2538, 24: 3040}
N26 = 20448  # 6*T_6 + 7*c(7,8)
N28 = 33240  # 7*c(7,7) + 8*c(8,10) + 9*c(9,13)
EXPECTED_TREES = {1: 1, 2: 3, 3: 15, 4: 83, 5: 486, 6: 2967, 7: 18748, 8: 121725, 9: 807381}
WAVE6_K0 = Fraction(848449, 500000)
LAMBDA = Fraction(11**11, 10**10)

DIRS = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))

# The 11-cell cavity witness: the six spokes around the origin plus five
# corner connectors spanning the spokes as the path
#   +x -(+y) -(-x) -(-y) -(+z)  with  (-y) -(-z),
# each connector a two-spoke corner cell (exactly two +-1 coordinates).
CAVITY_WITNESS = [
    (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1),
    (1, 1, 0), (-1, 1, 0), (-1, -1, 0), (0, -1, 1), (0, -1, -1),
]
# A 3x3x1 plate with the centre punched out: a tunnel, not a cavity.
TUNNEL_WITNESS = [(x, y, 0) for x in range(3) for y in range(3) if (x, y) != (1, 1)]
SPOKES = tuple(DIRS)


class StageBudgetExceeded(RuntimeError):
    pass


def cpu_time() -> float:
    return time.process_time()


def peak_rss_gib() -> float:
    raw = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    bytes_ = raw if sys.platform == "darwin" else raw * 1024
    return bytes_ / (1024**3)


def guard_rss() -> None:
    rss = peak_rss_gib()
    if rss > RSS_LIMIT_GIB:
        raise MemoryError(f"peak RSS {rss:.2f} GiB exceeds limit {RSS_LIMIT_GIB} GiB")


# ---------------------------------------------------------------------------
# Exact census engine (translations quotiented, orientations distinct)
# ---------------------------------------------------------------------------

def canonical_bytes(cells) -> bytes:
    """Sort the min-translated cells and pack each as three 4-bit coordinates."""
    mins = (min(c[0] for c in cells), min(c[1] for c in cells), min(c[2] for c in cells))
    ks = sorted(
        ((c[0] - mins[0]) << 8) | ((c[1] - mins[1]) << 4) | (c[2] - mins[2])
        for c in cells
    )
    out = bytearray()
    for k in ks:
        out += bytes((k >> 8, (k >> 4) & 15, k & 15))
    return bytes(out)


def canonical_key(cells) -> int:
    """Translation-min form as one Python integer.

    Each of the n sorted 12-bit min-translated keys becomes a 13-bit digit
    with the marker bit 0x1000 set, so the digit count n is recoverable from
    the integer alone (every digit is nonzero) and the encoding is injective.
    """
    mins = (min(c[0] for c in cells), min(c[1] for c in cells), min(c[2] for c in cells))
    ks = sorted(
        ((c[0] - mins[0]) << 8) | ((c[1] - mins[1]) << 4) | (c[2] - mins[2])
        for c in cells
    )
    acc = 0
    for k in ks:
        acc = (acc << 13) | 0x1000 | k
    return acc


def key_cells(data):
    """Decode a canonical form (int or bytes) back to cell triples."""
    if isinstance(data, int):
        out: list[tuple[int, int, int]] = []
        while data:
            d = data & 0x1FFF
            k = d & 0xFFF
            out.append((k >> 8, (k >> 4) & 15, k & 15))
            data >>= 13
        return out
    return [(data[i], data[i + 1], data[i + 2]) for i in range(0, len(data), 3)]


def grow_level(level: set) -> set:
    """All canonical (n+1)-cell sets obtained by adjoining one face neighbour."""
    out: set = set()
    add = out.add
    for data in level:
        cells = key_cells(data)
        members = set()
        for x, y, z in cells:
            members.add((x << 8) | (y << 4) | z)
        for (x, y, z) in cells:
            for dx, dy, dz in DIRS:
                nx, ny, nz = x + dx, y + dy, z + dz
                # the parent is normalised: only the shifted axis can go negative
                if nx < 0:
                    shift, newk = 1 << 8, (ny << 4) | nz
                elif ny < 0:
                    shift, newk = 1 << 4, (nx << 8) | nz
                elif nz < 0:
                    shift, newk = 1, (nx << 8) | (ny << 4)
                else:
                    newk = (nx << 8) | (ny << 4) | nz
                    if newk in members:
                        continue
                    shift = 0
                ks = [k + shift for k in members]
                ks.append(newk)
                ks.sort()
                if isinstance(data, int):
                    acc = 0
                    for k in ks:
                        acc = (acc << 13) | 0x1000 | k
                    add(acc)
                else:
                    buf = bytearray(len(ks) * 3)
                    j = 0
                    for k in ks:
                        buf[j] = k >> 8
                        buf[j + 1] = (k >> 4) & 15
                        buf[j + 2] = k & 15
                        j += 3
                    add(bytes(buf))
    return out


def shape_bonds(data) -> int:
    s = set()
    for x, y, z in key_cells(data):
        s.add((x << 8) | (y << 4) | z)
    b = 0
    for k in s:
        if (k + 256) in s:
            b += 1
        if (k + 16) in s:
            b += 1
        if (k + 1) in s:
            b += 1
    return b


def complement_connected(data) -> bool:
    """Exact flood fill in the one-cell-expanded bounding box."""
    cells = key_cells(data)
    xs = [c[0] for c in cells]
    ys = [c[1] for c in cells]
    zs = [c[2] for c in cells]
    x0, x1 = min(xs) - 1, max(xs) + 1
    y0, y1 = min(ys) - 1, max(ys) + 1
    z0, z1 = min(zs) - 1, max(zs) + 1
    Sx, W, D = x1 - x0 + 1, y1 - y0 + 1, z1 - z0 + 1
    YD = W * D
    solid = bytearray(Sx * YD)
    for (x, y, z) in cells:
        solid[((x - x0) * W + (y - y0)) * D + (z - z0)] = 1
    total = len(solid)
    stack = [0]
    solid[0] = 1
    cnt = 1
    while stack:
        p = stack.pop()
        x, r = divmod(p, YD)
        y, z = divmod(r, D)
        if x > 0:
            q = p - YD
            if not solid[q]:
                solid[q] = 1
                stack.append(q)
                cnt += 1
        if x < Sx - 1:
            q = p + YD
            if not solid[q]:
                solid[q] = 1
                stack.append(q)
                cnt += 1
        if y > 0:
            q = p - D
            if not solid[q]:
                solid[q] = 1
                stack.append(q)
                cnt += 1
        if y < W - 1:
            q = p + D
            if not solid[q]:
                solid[q] = 1
                stack.append(q)
                cnt += 1
        if z > 0:
            q = p - 1
            if not solid[q]:
                solid[q] = 1
                stack.append(q)
                cnt += 1
        if z < D - 1:
            q = p + 1
            if not solid[q]:
                solid[q] = 1
                stack.append(q)
                cnt += 1
    return cnt + len(cells) == total


def enumerate_census(n_max: int, budget_s: float):
    """BFS census on integer canonical keys with an exact cavity filter.

    For every translation class the complement is flood-filled in the
    one-cell-expanded bounding box: a class is a valid outer contour IFF its
    complement is connected (cavity-free).  Both strata (all classes and
    cavity-free classes) are tabulated; the certificate uses the latter.
    """
    t0 = cpu_time()
    level = {canonical_key([(0, 0, 0)])}
    level_counts = {1: 1}
    strata_all: dict[tuple[int, int], int] = {(1, 0): 1}
    strata_ok: dict[tuple[int, int], int] = {(1, 0): 1}
    n = 1
    while n < n_max:
        n += 1
        level = grow_level(level)
        level_counts[n] = len(level)
        level_cav = 0
        for data in level:
            b = shape_bonds(data)
            key = (n, b)
            strata_all[key] = strata_all.get(key, 0) + 1
            if complement_connected(data):
                strata_ok[key] = strata_ok.get(key, 0) + 1
            else:
                level_cav += 1
        print(
            f"  n={n}: {level_counts[n]} translation classes, {level_cav} with cavities "
            f"(cpu {cpu_time() - t0:.1f}s, rss {peak_rss_gib():.2f} GiB)",
            flush=True,
        )
        guard_rss()
        if cpu_time() - t0 > budget_s:
            raise StageBudgetExceeded(f"census exceeded {budget_s}s at n={n}")
    return level_counts, strata_all, strata_ok


def tree_census(n_max: int, budget_s: float) -> dict[int, int]:
    """Independent enumeration of tree polycubes by leaf extension.

    A tree polycube (adjacency graph a tree) with at least two cells has a
    leaf cell; conversely, adjoining a cell adjacent to exactly one existing
    cell preserves tree-ness, and every (n+1)-cell tree arises this way from
    some n-cell tree (delete a leaf).  Dedup by the same canonical form.
    """
    t0 = cpu_time()
    level = {canonical_key([(0, 0, 0)])}
    trees = {1: 1}
    for n in range(1, n_max):
        out: set[bytes] = set()
        add = out.add
        for data in level:
            cells = key_cells(data)
            members = set(cells)
            extras = list(cells)
            for (x, y, z) in cells:
                for dx, dy, dz in DIRS:
                    nx, ny, nz = x + dx, y + dy, z + dz
                    if (nx, ny, nz) in members:
                        continue
                    neighbours = sum(
                        (nx + d2[0], ny + d2[1], nz + d2[2]) in members for d2 in DIRS
                    )
                    # leaf extension: EXACTLY one bond to the old tree
                    if neighbours == 1:
                        add(canonical_bytes(extras + [(nx, ny, nz)]))
        level = out
        trees[n + 1] = len(level)
        print(
            f"  tree n={n + 1}: {len(level)} classes (cpu {cpu_time() - t0:.1f}s)",
            flush=True,
        )
        guard_rss()
        if cpu_time() - t0 > budget_s:
            raise StageBudgetExceeded(f"tree census exceeded {budget_s}s at n={n + 1}")
    return trees


def cells_connected(cells) -> bool:
    """Face-connectivity of a finite cell set by explicit BFS."""
    cells = set(cells)
    start = next(iter(cells))
    seen = {start}
    stack = [start]
    while stack:
        x, y, z = stack.pop()
        for dx, dy, dz in DIRS:
            c = (x + dx, y + dy, z + dz)
            if c in cells and c not in seen:
                seen.add(c)
                stack.append(c)
    return len(seen) == len(cells)


def cavity_frontier_search(max_extra: int = 5, budget_s: float = 900.0) -> dict[str, object]:
    """Exhaustive frontier enumeration of supersets of the 6 spokes.

    Any polycube P whose complement has a bounded component D with |D| = 1
    contains all six face-neighbours (spokes) of the enclosed cell.  Every
    CONNECTED P with spokes subset P and |P| <= 6 + max_extra arises by adding
    cells one at a time, each adjacent to the current set (reverse spanning-
    tree leaf order), so the frontier growth below generates all of them; the
    connectivity and cavity tests then filter exactly the polycubes.  With
    max_extra = 4 this certifies that no polycube with <= 10 cells has a
    single-cell cavity, and with max_extra = 5 it exhibits the minimum 11.
    """
    t0 = cpu_time()
    start = frozenset(SPOKES)
    current = {start}
    sizes: dict[int, int] = {0: 1}
    connected_sizes: dict[int, int] = {}
    min_cavity = None
    for step in range(1, max_extra + 1):
        nxt: set[frozenset] = set()
        for p in current:
            frontier = set()
            for (x, y, z) in p:
                for dx, dy, dz in DIRS:
                    c = (x + dx, y + dy, z + dz)
                    if c != (0, 0, 0) and c not in p:
                        frontier.add(c)
            for c in frontier:
                nxt.add(p | {c})
        current = nxt
        conn = [p for p in current if cells_connected(p)]
        sizes[step] = len(current)
        connected_sizes[step] = len(conn)
        if min_cavity is None:
            for p in conn:
                if not complement_connected(canonical_bytes(sorted(p))):
                    min_cavity = 6 + step
                    break
        print(
            f"  cavity frontier +{step}: {len(current)} supersets, "
            f"{len(conn)} connected (cpu {cpu_time() - t0:.1f}s)",
            flush=True,
        )
        if cpu_time() - t0 > budget_s:
            raise StageBudgetExceeded(f"cavity frontier exceeded {budget_s}s")
    return {
        "supersets_by_extra_cells": sizes,
        "connected_supersets_by_extra_cells": connected_sizes,
        "min_connected_cavity_cell_count": min_cavity if min_cavity is not None else -1,
        "searched_extra_cells": max_extra,
    }


# ---------------------------------------------------------------------------
# Exact rational enclosures of e^{-t} (no floating point anywhere)
# ---------------------------------------------------------------------------

def exp_neg_bounds(t: Fraction, order: int = EXP_ORDER) -> tuple[Fraction, Fraction]:
    """Exact rational enclosure lo <= e^{-t} <= hi for rational t >= 0.

    Write t = N*u with integer N >= 1 and u <= 1/2.  For 0 <= u <= 1/2 < 1
    the alternating Taylor series of e^{-u} has strictly decreasing terms, so
    every odd partial sum lies below e^{-u} and every even partial sum above;
    raising the bracket to the N-th power (monotone on positives) brackets
    e^{-t}.  All arithmetic is exact Fraction arithmetic.
    """
    if t < 0:
        raise ValueError("t must be >= 0")
    if t == 0:
        return Fraction(1), Fraction(1)
    n_sub = max(1, -(-2 * t.numerator // t.denominator))  # ceil(2t) >= 1
    u = t / n_sub
    if u > Fraction(1, 2):
        n_sub += 1
        u = t / n_sub
    term = Fraction(1)
    partial = Fraction(1)  # running partial sum S_j = sum_{i<=j} (-1)^i u^i/i!
    s_lo = partial  # S_1 will overwrite; S_0 = 1 is the coarsest upper bound
    s_hi = partial
    for j in range(1, order + 1):
        term = term * u / j
        partial = partial - term if j % 2 else partial + term
        if j % 2:
            s_lo = partial  # odd partial sums increase to e^{-u} from below
        else:
            s_hi = partial  # even partial sums decrease to e^{-u} from above
    if not (0 < s_lo < s_hi):
        raise AssertionError("alternating bracket malformed")
    return s_lo**n_sub, s_hi**n_sub


TRUNC_DEN = 10**48  # certified outward truncation grid for series evaluations


def outer_truncate(lo: Fraction, hi: Fraction, q: int = TRUNC_DEN) -> tuple[Fraction, Fraction]:
    """Round lo DOWN and hi UP to the grid 1/q: the true value stays inside.

    Keeps every downstream denominator ~q^A with bounded numerators, so the
    exact-Fraction series evaluations in the bisections below cost O(1) per
    term instead of paying gcd on ever-growing numerators; certificate
    direction is preserved because the rounding is strictly outward.
    """
    lo_scaled = lo * q
    hi_scaled = hi * q
    return (
        Fraction(lo_scaled.numerator // lo_scaled.denominator, q),
        Fraction(-((-hi_scaled.numerator) // hi_scaled.denominator), q),
    )


def series_T_coeffs(degree: int) -> list[int]:
    """Coefficients of T(z) solving T = z(1+T)^11 by iterative substitution.

    Starting from T = z, each pass T <- z(1+T)^11 stabilises one further
    coefficient, so degree+2 passes give the exact series through z^degree.
    """
    T = [0] * (degree + 1)  # T[i] = [z^i] T
    T[1] = 1
    for _ in range(degree + 2):
        base = [1] + T[1:]  # coefficients of 1 + T, index = power of z
        acc = [1] + [0] * degree
        for _power in range(11):
            nxt = [0] * (degree + 1)
            for i, a in enumerate(acc):
                if a:
                    for j, b in enumerate(base):
                        if b and i + j <= degree:
                            nxt[i + j] += a * b
            acc = nxt
        shifted = [0] * (degree + 1)
        for i in range(degree):
            shifted[i + 1] = acc[i]  # multiply by z
        if shifted == T:
            break
        T = shifted
    return T


def binomial(n: int, k: int) -> int:
    if k < 0 or k > n:
        return 0
    r = 1
    for i in range(k):
        r = r * (n - i) // (i + 1)
    return r


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def decimal_str(value: Fraction, places: int = 40) -> str:
    scale = 10**places
    scaled = value * scale
    lo = scaled.numerator // scaled.denominator
    hi = -((-scaled.numerator) // scaled.denominator)
    sign = "-" if lo < 0 else ""

    def fmt(integer: int) -> str:
        digits = str(abs(integer)).rjust(places + 1, "0")
        return f"{sign}{digits[:-places]}.{digits[-places:]}"

    return f"[{fmt(lo)}, {fmt(hi)}]"


def n_table(strata: dict[tuple[int, int], int]) -> dict[int, int]:
    """Exact N(A) for every area A reachable by cells n <= N_MAX."""
    table: dict[int, int] = {}
    for (n, b), count in strata.items():
        A = 6 * n - 2 * b
        table[A] = table.get(A, 0) + n * count
    return {A: c for A, c in table.items() if c > 0}


def main() -> None:
    t0 = cpu_time()
    checks: list[dict[str, object]] = []

    def check(name: str, passed: bool, detail: str) -> None:
        checks.append({"name": name, "passed": bool(passed), "detail": detail})
        print(f"  [{'ok' if passed else 'FAIL'}] {name}: {detail}", flush=True)

    # =====================================================================
    # Stage A: exact census through n = 10 cells
    # =====================================================================
    print(f"Stage A: census BFS through n={N_MAX} (complete for areas <= 28)", flush=True)
    ta = cpu_time()
    level_counts, strata_all, strata_ok = enumerate_census(N_MAX, BUDGET_CENSUS_S)
    stage_a_s = cpu_time() - ta
    guard_rss()
    for n in range(1, N_MAX + 1):
        total = sum(v for (nn, _), v in strata_all.items() if nn == n)
        check(f"level_sum_{n}", total == level_counts[n], f"{total} == {level_counts[n]}")
    for n, expected in EXPECTED_LEVELS.items():
        check(f"level_{n}", level_counts[n] == expected, f"{level_counts[n]} == {expected}")

    # flood-fill validity micro-tests
    check(
        "floodfill_cube",
        complement_connected(
            canonical_bytes([(x, y, z) for x in (0, 1) for y in (0, 1) for z in (0, 1)])
        ),
        "2x2x2 solid cube has connected complement",
    )
    check(
        "floodfill_tunnel",
        complement_connected(canonical_bytes(TUNNEL_WITNESS)),
        "3x3x1 plate with centre hole (tunnel) has connected complement",
    )
    cav_bonds = shape_bonds(canonical_bytes(CAVITY_WITNESS))
    check(
        "floodfill_cavity",
        (not complement_connected(canonical_bytes(CAVITY_WITNESS)))
        and cav_bonds == len(CAVITY_WITNESS) - 1,
        f"11-cell cavity tree detected (b={cav_bonds}=n-1, complement disconnected)",
    )
    n_cavities = sum(strata_all[k] - strata_ok.get(k, 0) for k in strata_all)

    # independent min-cavity frontier enumeration (|D| = 1 case) and the
    # |D| >= 2 neighbourhood lower bounds
    print("Stage B: min-cavity frontier enumeration around the six spokes", flush=True)
    tb = cpu_time()
    cavity = cavity_frontier_search(5, BUDGET_CAVITY_S)
    stage_b_s = cpu_time() - tb
    check(
        "min_cavity_is_11",
        cavity["min_connected_cavity_cell_count"] == 11
        and not complement_connected(canonical_bytes(CAVITY_WITNESS)),
        "no connected spoke-superset with <= 10 cells has a cavity and an "
 "11-cell connected one does: smallest single-cell-cavity polycube has "
        "exactly 11 cells",
    )
    # |D| = 2: P >= N(domino) (10 cells) and P = N(domino) is disconnected
    domino_nbhd = {(x + dx, y + dy, z + dz) for (x, y, z) in ((0, 0, 0), (1, 0, 0))
                   for (dx, dy, dz) in DIRS} - {(0, 0, 0), (1, 0, 0)}
    check(
        "domino_neighbourhood",
        len(domino_nbhd) == 10 and not cells_connected(domino_nbhd),
        f"|N(domino)| = {len(domino_nbhd)} and N(domino) is disconnected, so a "
        "two-cell cavity needs |P| >= 11",
    )
    # |D| >= 3: N(D) contains N(D3) for a 3-cell connected subset; minimum
    # over all 15 fixed three-cell classes (exhaustive)
    tri_min = min(
        len({(x + dx, y + dy, z + dz) for (x, y, z) in key_cells(data)
             for (dx, dy, dz) in DIRS} - set(key_cells(data)))
        for data in grow_level(grow_level({canonical_bytes([(0, 0, 0)])}))
    )
    check(
        "trio_neighbourhood_min",
        tri_min >= 11,
        f"min |N(D3)| over the 15 three-cell classes = {tri_min} >= 11, so every "
        "cavity with 3 <= |D| <= 5 needs |P| >= 16 - |D| >= 11 (larger |D| is "
        "excluded by the direct census computation below)",
    )
    check(
        "no_cavities_n_le_10",
        cavity["min_connected_cavity_cell_count"] == 11 and n_cavities == 0,
        "certificate: the exact flood fill ran over EVERY level-n class (n <= 10) "
        "and found 0 cavity-bearing classes (strata_all == strata_ok); the "
        "frontier search independently sharpens this to min cavity size = 11",
    )

    # wave-6 reproduction
    Ntab = n_table(strata_ok)
    check(
        "wave6_head",
        all(Ntab.get(A, 0) == WAVE6_HEAD[A] for A in WAVE6_HEAD),
        f"N(A), A<=24 = {[Ntab.get(A, 0) for A in sorted(WAVE6_HEAD)]}",
    )
    wave6_json = json.loads(WAVE6_PATH.read_text(encoding="utf-8"))
    stored_wave6 = {int(k): v for k, v in wave6_json["data"]["exact_head"]["N_of_A"].items()}
    check(
        "wave6_artifact_verbatim",
        all(Ntab.get(A, 0) == stored_wave6[A] for A in stored_wave6),
        "matches results/bounds/peierls_upper.json exact_head.N_of_A verbatim",
    )

    # Loomis-Whitney completeness cutoffs
    lw28 = 216 * 10**2 <= 28**3 < 216 * 11**2
    check(
        "loomis_whitney_28",
        lw28,
        "216*10^2=21600 <= 28^3=21952 < 26136=216*11^2: n<=10 complete for A<=28",
    )
    lw26 = 216 * 9**2 <= 26**3 < 216 * 10**2
    check(
        "loomis_whitney_26",
        lw26,
        "216*9^2=17496 <= 26^3=17576 < 21600: n<=9 complete for A<=26",
    )

    # area-26 and area-28 decompositions
    trees = {n: strata_ok[(n, n - 1)] for n in range(1, N_MAX + 1)}
    n26 = 6 * trees[6] + 7 * strata_ok[(7, 8)]
    check("N26", n26 == N26, f"6*{trees[6]} + 7*{strata_ok[(7, 8)]} = {n26}")
    n28 = 7 * strata_ok[(7, 7)] + 8 * strata_ok[(8, 10)] + 9 * strata_ok[(9, 13)]
    check(
        "N28",
        n28 == N28 and strata_ok.get((10, 16), 0) == 0,
        f"7*{strata_ok[(7, 7)]} + 8*{strata_ok[(8, 10)]} + 9*{strata_ok[(9, 13)]} = {n28}",
    )

    # tree cross-check by leaf extension
    print("Stage C: independent tree census by leaf extension", flush=True)
    tc = cpu_time()
    trees2 = tree_census(N_MAX, BUDGET_TREE_S)
    stage_c_s = cpu_time() - tc
    check(
        "tree_census",
        all(trees2[n] == trees[n] for n in range(1, N_MAX + 1)),
        f"T_n(leaf extension) == T_n(census stratum), n<=10: "
        f"{[trees2[n] for n in range(1, N_MAX + 1)]}",
    )
    check(
        "tree_expected_small",
        all(trees[n] == EXPECTED_TREES[n] for n in EXPECTED_TREES),
        "T_1..T_9 = 1, 3, 15, 83, 486, 2967, 18748, 121725, 807381",
    )
    check(
        "tree_closed_form",
        all(6 * n - 2 * (n - 1) == 4 * n + 2 for n in range(1, N_MAX + 1))
        and all(strata_all.get((n, b), 0) == 0 for n in range(1, N_MAX + 1) for b in range(n - 1)),
        "every tree class has area exactly A=4n+2=6n-2(n-1); no class has b<n-1",
    )

    # tail components M_A (n <= 10 parts of areas 28..42)
    parts: dict[int, dict[str, int]] = {}
    M: dict[int, int] = {}
    for A in TAIL_AREAS:
        comps = {}
        tot = 0
        for n in range(1, N_MAX + 1):
            b = 3 * n - A // 2
            c = strata_ok.get((n, b), 0)
            if c:
                comps[str(n)] = c
                tot += n * c
        M[A] = tot
        parts[A] = comps
    check(
        "tail_components",
        M[28] == N28
        and M[30] == 7 * trees[7] + 8 * strata_ok[(8, 9)] + 9 * strata_ok[(9, 12)] + 10 * strata_ok[(10, 15)]
        and M[42] == 10 * trees[10],
        f"M_A (n<=10 parts) = {[M[A] for A in TAIL_AREAS]}",
    )

    # =====================================================================
    # Stage D: exact rational threshold / residual arithmetic
    # =====================================================================
    print("Stage D: exact rational enclosures and certificates", flush=True)
    td = cpu_time()
    from mpmath import iv as _iv
    import mpmath as _mp

    _iv.dps = IV_DPS
    I3_iv = (
        _iv.sqrt(6)
        * _iv.gamma(_iv.mpf(1) / 24)
        * _iv.gamma(_iv.mpf(5) / 24)
        * _iv.gamma(_iv.mpf(7) / 24)
        * _iv.gamma(_iv.mpf(11) / 24)
        / (96 * _iv.pi**3)
    )
    with _mp.workdps(2 * IV_DPS):

        def iv_frac(v) -> Fraction:
            f = _mp.mpf(v)
            sign, man, exp, _bc = f._mpf_
            if man == 0:
                return Fraction(0)
            out = Fraction(man) * Fraction(2) ** exp
            return -out if sign else out

        i3_iv_lo = iv_frac(I3_iv.a)
        i3_iv_hi = iv_frac(I3_iv.b)

    # repo-certified anchor interval (proofs/upper_infrared.md chain)
    from decimal import Decimal

    infrared = json.loads(INFRARED_PATH.read_text(encoding="utf-8"))
    i3_repo_lo = Fraction(Decimal(infrared["data"]["certified_constants"]["I3"][0]))
    i3_repo_hi = Fraction(Decimal(infrared["data"]["certified_constants"]["I3"][1]))
    eps = Fraction(1, 10**90)
    check(
        "I3_anchor_nesting",
        i3_repo_lo - eps <= i3_iv_lo and i3_iv_hi <= i3_repo_hi + eps,
        "fresh mpmath dps=100 I_3 interval nested in the repo-certified interval "
        "of results/bounds/upper_infrared.json (widened by 1e-90)",
    )
    # Proof-side enclosure of I_3: the repo-certified interval itself.
    i3_lo = max(i3_iv_lo, i3_repo_lo)
    i3_hi = min(i3_iv_hi, i3_repo_hi)
    check("I3_intersection", i3_lo < i3_hi, f"I_3 in {decimal_str(i3_lo, 40)}")
    # x* = e^{-I_3} with I_3 in [i3_lo, i3_hi]: monotonicity gives the bracket
    e1_lo, e1_hi = exp_neg_bounds(i3_lo)  # encloses e^{-i3_lo}  (upper end)
    e2_lo, e2_hi = exp_neg_bounds(i3_hi)  # encloses e^{-i3_hi}  (lower end)
    x_lo, x_hi = e2_lo, e1_hi
    with _mp.workdps(2 * IV_DPS):
        x_iv = _iv.exp(-I3_iv)
        x_iv_lo = iv_frac(x_iv.a)
        x_iv_hi = iv_frac(x_iv.b)
    check(
        "x_interval_agreement",
        x_lo < x_hi and x_iv_lo <= x_hi and x_lo <= x_iv_hi,
        f"pure-Fraction and mpmath x*=e^{{-I_3}} enclosures overlap: "
        f"{decimal_str(x_lo, 30)}",
    )
    check(
        "x_interval_width",
        (x_hi - x_lo) < Fraction(1, 10**60),
        f"width < 1e-60 ({decimal_str(x_hi - x_lo, 12)})",
    )

    head = {A: Ntab[A] for A in (6, 10, 14, 16, 18, 20, 22, 24, 26)}
    head28 = dict(head)
    head28[28] = Ntab[28]

    def series_value(coeffs: dict[int, int], y: Fraction) -> Fraction:
        # ascending powers: one multiplication per term instead of one y**A
        total = Fraction(0)
        power = Fraction(1)
        a_prev = 0
        for A in sorted(coeffs):
            power *= y ** (A - a_prev)
            a_prev = A
            total += coeffs[A] * power
        return total

    # ---- master limitation inequality at K_bar ----
    H26_lo = series_value(head, x_lo)
    H26_hi = series_value(head, x_hi)
    H28_lo = series_value(head28, x_lo)
    H28_hi = series_value(head28, x_hi)
    tail_lower = series_value(M, x_lo)  # <= true tail mass at areas >= 28
    S_lo = H26_lo + tail_lower  # <= V(x*)
    check(
        "master_inequality",
        S_lo > Fraction(1, 2),
        f"V(x*) >= head26(x_lo)+tail(28..42)(x_lo) > 1/2; certified margin "
        f"{decimal_str(S_lo - Fraction(1, 2), 20)}",
    )
    # same-endpoint comparison: tail(x*) >= tail(x_lo) and slack(x*) <= 1/2 - H26(x_lo)
    slack_at_lo = Fraction(1, 2) - H26_lo
    check(
        "tail_exceeds_slack",
        tail_lower > slack_at_lo,
        f"certified tail mass at areas>=28 {decimal_str(tail_lower, 12)} > certified "
        f"slack upper bound 1/2-head26(x_lo) {decimal_str(slack_at_lo, 12)}",
    )
    V28_lo = H28_lo + (tail_lower - Fraction(M[28]) * x_lo**28)

    # ---- K* bisection: threshold of the certified lower bound V ----
    def V_lower_at_K(K: Fraction) -> Fraction:
        y_lo, _ = outer_truncate(*exp_neg_bounds(2 * K))
        return series_value(head28, y_lo) + (
            series_value(M, y_lo) - Fraction(M[28]) * y_lo**28
        )

    def V_upper_at_K(K: Fraction) -> Fraction:
        _, y_hi = outer_truncate(*exp_neg_bounds(2 * K))
        return series_value(head28, y_hi) + (
            series_value(M, y_hi) - Fraction(M[28]) * y_hi**28
        )

    kbar_lo, kbar_hi = i3_lo / 2, i3_hi / 2
    check(
        "V_at_kbar",
        V_lower_at_K(kbar_hi) > Fraction(1, 2),
        "V(e^{-2K}) > 1/2 already certified at K = upper end of the I_3/2 interval",
    )
    k_lo, k_hi = kbar_hi, Fraction(3, 2)  # V(k_lo side) > 1/2 > V(k_hi side)
    assert V_lower_at_K(k_lo) > Fraction(1, 2) and V_upper_at_K(k_hi) < Fraction(1, 2)
    for _ in range(72):
        mid = (k_lo + k_hi) / 2
        if V_lower_at_K(mid) > Fraction(1, 2):
            k_lo = mid  # certified V(mid) > 1/2: threshold lies above mid
        elif V_upper_at_K(mid) < Fraction(1, 2):
            k_hi = mid  # certified V(mid) < 1/2: threshold lies below mid
        else:
            break  # certified endpoints straddle 1/2: bracket already tight
    check(
        "K_star_bracket",
        k_hi - k_lo < Fraction(1, 10**18)
        and V_lower_at_K(k_lo) > Fraction(1, 2)
        and V_upper_at_K(k_hi) < Fraction(1, 2),
        f"K* in {decimal_str(k_lo, 18)}; V(e^{{-2K}}) > 1/2 for every K <= K*_lo, "
        "so the plain certificate fails for every K <= K*",
    )
    gap_lo = k_lo - kbar_hi  # > 0
    check(
        "K_star_above_incumbent",
        gap_lo > 0,
        f"certified limitation gap K*_lo - K_bar_hi = {decimal_str(gap_lo, 12)} > 0",
    )

    # ---- spurious head-only crossing K_26 and its refutation ----
    def H26_upper_at_K(K: Fraction) -> Fraction:
        _, y_hi = outer_truncate(*exp_neg_bounds(2 * K))
        return series_value(head, y_hi)

    def H26_lower_at_K(K: Fraction) -> Fraction:
        y_lo, _ = outer_truncate(*exp_neg_bounds(2 * K))
        return series_value(head, y_lo)

    b_lo, b_hi = Fraction(1, 8), kbar_lo
    assert H26_upper_at_K(b_lo) > Fraction(1, 2) and H26_lower_at_K(b_hi) < Fraction(1, 2)
    for _ in range(72):
        mid = (b_lo + b_hi) / 2
        if H26_upper_at_K(mid) < Fraction(1, 2):
            b_hi = mid  # certified head < 1/2 at mid
        elif H26_lower_at_K(mid) > Fraction(1, 2):
            b_lo = mid  # certified head > 1/2 at mid
        else:
            break  # enclosure straddles 1/2: certified bracket already tight
    check(
        "K26_spurious_bracket",
        H26_upper_at_K(b_lo) > Fraction(1, 2) and H26_lower_at_K(b_hi) < Fraction(1, 2),
        f"head-only crossing K_26 in {decimal_str(b_lo, 18)} (deep below K_bar)",
    )
    refuted_margin = V_lower_at_K(b_hi)
    check(
        "K26_refuted",
        refuted_margin > Fraction(1, 2),
        f"at the spurious crossing the full certified lower bound is still "
        f"{decimal_str(refuted_margin, 12)} > 1/2: the n=26 truncation feigns a "
        f"gain of K_bar - K_26 = {decimal_str(kbar_lo - b_hi, 12)} that the proved "
        "tail mass annihilates",
    )

    # ---- machine re-proof of the wave-6 exponential tail bound ----
    degree = 33
    Tc = series_T_coeffs(degree)
    # R = z (1+T)^12 coefficients
    R = [0] * (degree + 1)
    base = [1] + Tc[1:]  # coefficients of 1 + T, index = power of z
    acc = [1] + [0] * degree
    for _ in range(12):
        nxt = [0] * (degree + 1)
        for i, a in enumerate(acc):
            if a:
                for j, b in enumerate(base):
                    if b and i + j <= degree:
                        nxt[i + j] += a * b
        acc = nxt
    for i in range(degree):
        R[i + 1] = acc[i]
    lagrange_ok = R[1] == 1
    for A in range(2, degree + 1):
        lhs = R[A]
        rhs = 12 * binomial(11 * A, A - 2) // (A - 1)
        lagrange_ok = lagrange_ok and lhs == rhs
    check(
        "lagrange_identity",
        lagrange_ok,
        "[z^{A-1}] z(1+T)^12 == 12/(A-1) C(11A, A-2) for A=2..33 by iterative "
        "series reversion of T = z(1+T)^11",
    )
    binom_ok = True
    for A in range(26, 201, 2):
        binom_ok = binom_ok and 100 * binomial(11 * A, A - 2) <= LAMBDA**A
    check(
        "binomial_majorant_rows",
        binom_ok,
        "100*C(11A, A-2) <= lambda^A for every even A in [26, 200] (exact integers)",
    )
    animal_ok = True
    for A in range(26, 201, 2):
        B_A = 12 * binomial(11 * A, A - 2) // (A - 1)
        animal_ok = animal_ok and 625 * A * B_A <= 78 * LAMBDA**A
    check(
        "animal_bound_rows",
        animal_ok,
        "625*A*B_A <= 78*lambda^A for every even A in [26, 200], i.e. N(A) <= "
        "(78/625) lambda^A (root multiplicity <= A)",
    )

    # ---- certified head+lambda-tail grid endpoint ----
    def cert_upper_at_K(K: Fraction) -> Fraction:
        _, y = outer_truncate(*exp_neg_bounds(2 * K))  # y = upper end >= e^{-2K}
        ly = LAMBDA * y
        if ly >= 1:
            return Fraction(10**6)
        h = series_value(head28, y)
        tail = Fraction(78, 625) * ly**30 / (1 - ly * ly)
        return h + tail

    def cert_lower_at_K(K: Fraction) -> Fraction:
        y, _ = outer_truncate(*exp_neg_bounds(2 * K))
        ly = LAMBDA * y
        h = series_value(head28, y)
        tail = Fraction(78, 625) * ly**30 / (1 - ly * ly)
        return h + tail

    # bracket: certified above at 1.69, certified below at 7/4; the crossing
    # (~1.6954) lies inside, strictly below the wave-6 endpoint 1.696898
    c_lo, c_hi = Fraction(169, 100), Fraction(7, 4)
    assert cert_lower_at_K(c_lo) > Fraction(1, 2)
    assert cert_upper_at_K(c_hi) < Fraction(1, 2)
    assert cert_lower_at_K(WAVE6_K0) < Fraction(1, 2)  # wave-6 endpoint passes
    for _ in range(60):
        mid = (c_lo + c_hi) / 2
        if cert_lower_at_K(mid) > Fraction(1, 2):
            c_lo = mid  # certified above 1/2: crossing lies above mid
        elif cert_upper_at_K(mid) < Fraction(1, 2):
            c_hi = mid  # certified below 1/2: crossing lies below mid
        else:
            break  # enclosure straddles 1/2: bracket already certified
    grid = 10**6
    g = Fraction(-(-c_hi.numerator * grid // c_hi.denominator), grid)  # ceil to 1e-6
    g_pass = cert_upper_at_K(g) < Fraction(1, 2)
    g_below_fail = cert_lower_at_K(g - Fraction(1, grid)) > Fraction(1, 2)
    check(
        "lambda_tail_grid_endpoint",
        g_pass and g_below_fail and g < WAVE6_K0,
        f"least 1e-6-grid K with head28+(78/625)(ly)^30/(1-(ly)^2) < 1/2 is "
        f"{fraction_text(g)} = {decimal_str(g, 8)}; improves wave-6 K_0 = "
        f"{fraction_text(WAVE6_K0)} by {fraction_text(WAVE6_K0 - g)}",
    )

    # ---- convergence threshold of the exponential tail vs the incumbent ----
    def log_lambda_half_bisect() -> tuple[Fraction, Fraction]:
        lo, hi = Fraction(0), Fraction(4)
        for _ in range(72):
            mid = (lo + hi) / 2
            y_lo, y_hi = outer_truncate(*exp_neg_bounds(2 * mid))
            if LAMBDA * y_hi < 1:  # certified lambda e^{-2 mid} < 1
                hi = mid
            elif LAMBDA * y_lo > 1:  # certified lambda e^{-2 mid} > 1
                lo = mid
            else:
                break
        return lo, hi

    conv_lo, conv_hi = log_lambda_half_bisect()
    check(
        "tail_convergence_threshold",
        conv_lo > kbar_hi,
        f"(1/2) log lambda in {decimal_str(conv_lo, 12)}: the exponential tail "
        f"converges only above this, a certified {decimal_str(conv_lo - kbar_hi, 12)} "
        "above K_bar; no finite exact head can close this gap (wave-6 Section 5)",
    )

    stage_d_s = cpu_time() - td
    guard_rss()

    # =====================================================================
    # Artifact
    # =====================================================================
    strata_serial = {
        f"{n},{b}": {"all": strata_all[(n, b)], "cavity_free": strata_ok.get((n, b), 0)}
        for (n, b) in sorted(strata_all)
    }
    result = {
        "meta": {
            "script": SCRIPT,
            "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%S+00:00", time.gmtime()),
            "classification": (
                "[COMPUTATION] exact rooted outer-contour census, areas <= 28 complete "
                "(cells n <= 10, Loomis-Whitney cutoff), bond strata n <= 10; "
                "[THEOREM] the plain Peierls contour-sum certificate fails for every "
                "K <= K* with certified K* > I_3/2, so the exact-count head-extension "
                "program cannot reach or beat the incumbent upper endpoint I_3/2; "
                "[COMPUTATION] machine-reproved exponential tail gives the least "
                "1e-6-grid head+tail endpoint, improving wave-6's 1.696898"
            ),
            "arithmetic": {
                "combinatorics": "exact Python integers",
                "transcendental": (
                    "pure exact-Fraction alternating-Taylor enclosures of e^{-2K} "
                    "(subdivision u <= 1/2), cross-checked against mpmath interval "
                    "arithmetic at dps=100 (proofs/upper_infrared.md Section 2); "
                    "every comparison decided on exact Fraction endpoints"
                ),
                "floating_point_used_in_certificate": False,
            },
            "constraint_class": (
                "Peierls outer-contour counting class only; the audited "
                "reflection-positivity/infrared two-point class is exactly optimal "
                "at I_3/2 (proofs/upper_infrared.md) and its zero-magnetization "
                "profile saturates the constraint (proofs/mag_floor.md), so it is "
                "ruled out; no two-point-only variant is offered"
            ),
            "stage_budgets_cpu_s": {
                "census": BUDGET_CENSUS_S,
                "cavity_frontier": BUDGET_CAVITY_S,
                "tree_census": BUDGET_TREE_S,
                "arithmetic": BUDGET_ARITH_S,
            },
            "stage_cpu_seconds": {
                "census": round(stage_a_s, 1),
                "cavity_frontier": round(stage_b_s, 1),
                "tree_census": round(stage_c_s, 1),
                "arithmetic": round(stage_d_s, 1),
            },
            "peak_rss_gib": round(peak_rss_gib(), 3),
            "provenance": SCRIPT,
            "cpu_seconds": round(cpu_time() - t0, 1),
        },
        "data": {
            "model": "nearest-neighbour ferromagnetic Ising model on Z^3, K=beta*J, x=e^{-2K}",
            "contour_definition": (
                "edge-connected outer dual-plaquette boundary of the filled minus "
                "component containing the origin (proofs/kc_upper_peierls.md Sections "
                "2-3); N(A) = sum over translation classes with n cells of weight n"
            ),
            "census": {
                "cell_cutoff": N_MAX,
                "complete_for_areas_upto": 28,
                "level_counts": {str(n): level_counts[n] for n in sorted(level_counts)},
                "bond_strata": strata_serial,
                "all_strata_cavity_free": n_cavities == 0,
                "cavity_floodfill_microtests": [
                    "2x2x2 cube -> connected",
                    "3x3x1 plate with centre hole (tunnel) -> connected",
                    f"11-cell cavity tree (b={cav_bonds}) -> disconnected, detected",
                ],
                "min_cavity_search": cavity,
                "min_cavity_cell_count": cavity["min_connected_cavity_cell_count"],
                "cavity_free_certificate": (
                    "exact complement flood fill over every level class n<=10: "
                    "0 cavity-bearing classes"
                ),
            },
            "tree_classes": {
                "identity": "tree class: b = n-1 bonds, area A = 6n - 2(n-1) = 4n + 2 exactly",
                "T_n": {str(n): trees[n] for n in sorted(trees)},
                "T_n_leaf_extension_crosscheck": {str(n): trees2[n] for n in sorted(trees2)},
                "successive_ratios": {
                    str(n): fraction_text(Fraction(trees[n + 1], trees[n]))
                    for n in range(1, N_MAX)
                },
            },
            "exact_N_table": {str(A): Ntab[A] for A in sorted(Ntab) if A <= 28},
            "N26_decomposition": {"6*T_6": 6 * trees[6], "7*c(7,8)": 7 * strata_ok[(7, 8)]},
            "N28_decomposition": {
                "7*c(7,7)": 7 * strata_ok[(7, 7)],
                "8*c(8,10)": 8 * strata_ok[(8, 10)],
                "9*c(9,13)": 9 * strata_ok[(9, 13)],
                "10*c(10,16)": 0,
            },
            "tail_lower_bound_components": {
                str(A): {"M_A_n_le_10": M[A], "strata": parts[A]} for A in TAIL_AREAS
            },
            "limitation_certificate": {
                "kbar_interval": {
                    "lower_exact": fraction_text(kbar_lo),
                    "upper_exact": fraction_text(kbar_hi),
                    "lower_decimal": decimal_str(kbar_lo, 40),
                    "upper_decimal": decimal_str(kbar_hi, 40),
                },
                "x_star_interval": {
                    "lower_exact": fraction_text(x_lo),
                    "upper_exact": fraction_text(x_hi),
                    "lower_decimal": decimal_str(x_lo, 50),
                    "upper_decimal": decimal_str(x_hi, 50),
                },
                "head26_at_x_lo_decimal": decimal_str(H26_lo, 30),
                "head26_at_x_hi_decimal": decimal_str(H26_hi, 30),
                "head28_at_x_lo_decimal": decimal_str(H28_lo, 30),
                "head28_at_x_hi_decimal": decimal_str(H28_hi, 30),
                "tail_lower_decimal": decimal_str(tail_lower, 30),
                "master_sum_lower_decimal": decimal_str(S_lo, 30),
                "margin_over_one_half_decimal": decimal_str(S_lo - Fraction(1, 2), 30),
                "margin_positive_certified_on_exact_fractions": True,
                "slack_upper_bound_decimal": decimal_str(slack_at_lo, 30),
                "K_star_interval_exact": [fraction_text(k_lo), fraction_text(k_hi)],
                "K_star_decimal": decimal_str(k_lo, 20),
                "K_star_minus_kbar_exact": fraction_text(gap_lo),
                "K_star_minus_kbar_decimal": decimal_str(gap_lo, 20),
                "K26_spurious_interval_exact": [fraction_text(b_lo), fraction_text(b_hi)],
                "K26_spurious_decimal": decimal_str(b_lo, 20),
                "K26_refutation_value_decimal": decimal_str(refuted_margin, 20),
                "K26_feigned_gain_kbar_minus_K26_decimal": decimal_str(kbar_lo - b_hi, 20),
                "statement": (
                    "V(x) := sum_{A<=26} N(A) x^A + sum_{A in 28..42 even} M_A x^A is a "
                    "certified lower bound of the full contour sum at x = e^{-2K}; "
                    "V(x*) > 1/2 and V(e^{-2K}) > 1/2 for every K <= K* (K* > I_3/2 "
                    "certified), so the plain Peierls certificate of "
                    "kc_upper_peierls.md (2.5)-(2.6) fails for every K <= K*"
                ),
            },
            "proved_exponential_tail": {
                "lambda_exact": fraction_text(LAMBDA),
                "bound": "N(A) <= A*B_A, B_A = 12/(A-1) C(11A, A-2); 100*C(11A,A-2) "
                "<= lambda^A; hence N(A) <= (78/625) lambda^A for A >= 26",
                "verified_rows": "A = 2..33 (Lagrange identity, iterative reversion); "
                "even A in [26,200] (binomial majorant and animal bound, exact)",
                "tail_formula": "sum_{A>=30, even} N(A) y^A <= (78/625) (lambda y)^30 / (1-(lambda y)^2)",
                "grid_endpoint_exact": fraction_text(g),
                "grid_endpoint_decimal": decimal_str(g, 8),
                "grid_endpoint_minus_one_fails": True,
                "wave6_endpoint_exact": fraction_text(WAVE6_K0),
                "improvement_over_wave6_exact": fraction_text(WAVE6_K0 - g),
                "convergence_threshold_half_log_lambda": [
                    fraction_text(conv_lo),
                    fraction_text(conv_hi),
                ],
                "convergence_threshold_decimal": decimal_str(conv_lo, 16),
                "convergence_minus_kbar_decimal": decimal_str(conv_lo - kbar_hi, 16),
                "note": (
                    "an honest but weak improvement of the wave-6 certified endpoint; "
                    "still far above the incumbent I_3/2 = 0.2527, and the tail "
                    "converges only above (1/2) log lambda = 1.6754..., so no finite "
                    "head extension of this tail can approach the incumbent"
                ),
            },
            "open": (
                "whether the full contour series sum_A N(A) x*^A converges or diverges "
                "at x* = e^{-I_3} remains open; the tree stratum alone converges iff "
                "lim T_{n+1}/T_n < x*^{-4}; successive tree ratios through n=10 are "
                + ", ".join(fraction_text(Fraction(trees[n + 1], trees[n])) for n in range(6, N_MAX))
                + " (rising, last "
                + decimal_str(Fraction(trees[N_MAX], trees[N_MAX - 1]), 4)
                + ") versus x*^{-4} in "
                + decimal_str(Fraction(1) / x_hi**4, 6)
                + " .. "
                + decimal_str(Fraction(1) / x_lo**4, 6)
            ),
        },
        "checks": checks,
    }
    RESULT_PATH.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")

    failed = [c["name"] for c in checks if not c["passed"]]
    print(f"\nWrote {RESULT_PATH}", flush=True)
    print(f"cpu seconds: {cpu_time() - t0:.1f}, peak RSS {peak_rss_gib():.2f} GiB", flush=True)
    print(f"limitation margin over 1/2 at K_bar: {decimal_str(S_lo - Fraction(1, 2), 20)}", flush=True)
    print(f"K* (certificate fails below): {decimal_str(k_lo, 16)}", flush=True)
    print(f"K* - K_bar: {decimal_str(gap_lo, 16)}", flush=True)
    if failed:
        raise SystemExit(f"FAIL: {failed}")
    print("PASS: census + limitation certificate + tail endpoint regenerated", flush=True)


if __name__ == "__main__":
    main()
