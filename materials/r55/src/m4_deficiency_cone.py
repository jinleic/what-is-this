"""The m=4 high-order deficiency cone over R(4,5) strata.

Task 2 of notes/higher_order_identity_m4_plan_2026-08-17.md.  This module
extends the frozen m=3 cone (`m3_deficiency_cone`, whose domain and
g-intervals are reused verbatim by import) with the five streaming motif
aggregates of the landed m=4 kernel (`m4_subgraph_identities`):

    q       = sum_x C(d_x, 2) - t        (the m=3 wedge aggregate)
    t       = triangle count
    diamond = induced K4 minus an edge
    k3k1    = induced K3 + K1
    i4      = independent 4-sets

Every aggregate is a linear combination of the eleven induced 4-vertex motif
counts, so ALL of {q, t, diamond, k3k1, i4} are constrained by the rank
relations R1-R5 below.  The relaxed envelope LP over those relations yields
a certified outer window for any (order, edges) class whose complete catalog
is unavailable; the relaxed right-hand sides aggregate only degree histograms
whose per-neighborhood R(3,5) edge windows survive, exactly as the frozen m=3
relaxation does (same feasibility semantics, exact integer arithmetic).

Design decisions pinned by the plan and the TDD suite
(tests/test_m4_deficiency_cone.py):

* Envelope LP: {s >= 0 : R1, R2 exact; R3, R4, R5 two-sided}; its vertices
  are enumerated as feasible BASES over the 17 trivial/range rows; the
  exactness-preserving cut is that any basis activating both endpoints of
  the same range row is singular, and any basis with active > 3 must do so
  by pigeonhole - so only zeroed_count >= 6 strata reach square solves, and
  each with at most one endpoint per range row.  Counts and per-key min/max
  are provably identical to the naive C(17, 9) enumeration, and a zeroed
  pick coordinate is a legitimate vertex (offered as Fraction(0)).
* Envelope failures are fail-closed: ANY exception in the LP path records
  the class in `lp_fallback_classes()` and yields the trivially sound
  window (0, C(order, 4)).
* State h intervals use DESIGN A (frozen m=3 convention): the y-side i4
  window is resolved at the PLAIN stratum class (m, E45(m) - b) with
  m = N_TARGET - 1 - d, matching the m=3 q(D_v) convention.
* Self-loading policy: the g side reproduces the 3,215 frozen m=3 records
  exactly (7-field equality by construction); h motif windows resolve from
  data/higher_identity_m4.json when it exists, else the sound trivial
  [0, C(r, 4)] window.  Focused tests run in seconds either way.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import os
import sys
import tempfile
import time
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

from check_ramsey import parse_graph6_line, popcount
from subgraph_identities import triangle_count, complement, edge_count
from m3_deficiency_cone import (
    InfeasibleEdgeClass,
    degree_bounds,
    degree_histograms,
    load_edge_bounds,
    r35_edge_windows,
    _edge_window,
    _count_graph6_lines,
    build_states as _m3_build_states,
    m3_state_interval,
    _select_main_records,
    _select_extreme_records,
    _input_record,
    _verify_unchanged_input,
    _sorted_inputs,
    _relative_path,
    _load_json,
    _exact_int,
    R35_FILES,
    R45_CENSUS_FILE,
    R55_REPLAY_FILE,
    R55_REPLAY_GRAPHS,
    R55_REPLAY_ORDER,
    R45_EXTREME_DIRNAME,
    TRUST_ROOTS as M3_TRUST_ROOTS,
    STATE_DEGREES,
    N_TARGET,
    DATA_DIRNAME,
    TABLE_FILE,
    VALIDATION_FILE,
    VALIDATION_EXTREME_FILE,
    OUTER_SOURCE,
)
from m4_subgraph_identities import (
    diamond_count,
    induced_k3k1_count,
    independent_quad_count,
    m4_residual_specialized,
)

_SRC_DIR = Path(__file__).resolve().parent
_DATA_DIR = _SRC_DIR.parent / DATA_DIRNAME
_EXTREME_DIR = _DATA_DIR / "r45extreme"
_M3_ARTIFACT_PATH = _DATA_DIR / "higher_identity_m3.json"
_M4_ARTIFACT_PATH = _DATA_DIR / "higher_identity_m4.json"

SCHEMA_VERSION = 2
CAMPAIGN_ID = "higher_order_identity_positive_deficiency_m4"
DISPOSITION = "M4_ANALYZED"        # Task 3 owns every later disposition
MOTIF4_CATALOG_SOURCE = "catalog:higher_identity_m4.json"
TRIVIAL_SOURCE = "trivial:binomial_coefficient_bound"
TRUST_ROOTS = M3_TRUST_ROOTS + (
    (
        "The m=4 identity is locally re-derived (module docstring counting "
        "ledger of m4_subgraph_identities.py), exhaustively verified on all "
        "labeled graphs of order <= 6, spot-verified at orders 7-8, "
        "replayed on all 656 published n=42 graphs and complements, and "
        "replayed against the MR97 Theorem 3.2 n=49 numbers; primary-source "
        "cross-check anchors are MR97 Theorem 2.2 (m=4 case) and "
        "Theorem 3.2.",
        "https://users.cecs.anu.edu.au/~bdm/papers/r55.pdf",
    ),
    (
        "The local motif-window tables in data/higher_identity_m4.json are "
        "derived only from hash-pinned frozen catalogs plus the relation "
        "system R1-R5 (constants unit-tested) and are consumed downstream "
        "only through that artifact.",
        "https://arxiv.org/abs/2409.15709",
    ),
)

# The eleven induced 4-vertex motifs in the plan's canonical order.
MOTIF_ORDER = (
    "E4", "K2+2K1", "2K2", "P3+K1", "P4", "K1_3",
    "K3+K1", "C4", "paw", "diamond", "K4",
)
# Motif indices of the three raw-coordinate LP picks (diamond, K3+K1, E4).
_PICK_INDICES = (9, 6, 0)
MOTIF4_KEYS = ("q", "t", "diamond", "k3k1", "i4")

_MOTIF_EDGES = (0, 1, 2, 2, 3, 3, 3, 4, 4, 5, 6)
_MOTIF_DEGSQ = (0, 2, 4, 6, 10, 12, 12, 16, 18, 26, 36)
_MOTIF_TAU = (0, 0, 0, 0, 0, 0, 1, 0, 1, 2, 4)

# i4-index inside MOTIF4_KEYS is unused; motif windows are keyed by name.


def _c2(n):
    return n * (n - 1) // 2


def _c3(n):
    return n * (n - 1) * (n - 2) // 6


def _c4(n):
    return n * (n - 1) * (n - 2) * (n - 3) // 24


def motif4_vector(adj):
    """The five streaming aggregates of one graph, from the m=4 kernel."""
    triangles = triangle_count(adj)
    wedge = sum(_c2(popcount(row)) for row in adj)
    return {
        "q": wedge - triangles,
        "t": triangles,
        "diamond": diamond_count(adj),
        "k3k1": induced_k3k1_count(adj),
        "i4": independent_quad_count(adj),
    }


def affine_interval(coef, window):
    """Sign-aware scaling of one exact interval by an integer coefficient."""
    lo, hi = window
    if coef >= 0:
        return (coef * lo, coef * hi)
    return (coef * hi, coef * lo)


def relation_rows():
    """The R1-R5 coefficient rows over MOTIF_ORDER (plan constants verbatim)."""
    r2 = _MOTIF_EDGES
    return {
        "R1": (1,) * 11,
        "R2": r2,
        "R3": tuple(e * e for e in r2),
        "R4": _MOTIF_DEGSQ,
        "R5": _MOTIF_TAU,
    }


def _f_of(order, degree):
    """R4 right-side per-degree term f(d) = d*C(r-1-d,2)+4*C(d,2)(r-1-d)+9*C(d,3)."""
    return (degree * _c2(order - 1 - degree)
            + 4 * _c2(degree) * (order - 1 - degree)
            + 9 * _c3(degree))


def _solve_square_int(mat):
    """Exact solution of one square integer system, or None if singular.

    Fraction-free (Bareiss) forward elimination -- every intermediate entry
    stays an exact minor, so each division is exact -- followed by one
    Fraction back-substitution.  Mirrors the in-test reference semantics:
    unique-solution-or-None, so zero-count equality is solver-agnostic.
    """
    n = len(mat)
    a = [row[:] for row in mat]
    prev = 1
    for k in range(n - 1):
        if a[k][k] == 0:
            replacement = None
            for i in range(k + 1, n):
                if a[i][k]:
                    replacement = i
                    break
            if replacement is None:
                return None
            a[k], a[replacement] = a[replacement], a[k]
        pivot = a[k][k]
        top = a[k]
        for i in range(k + 1, n):
            row = a[i]
            factor = row[k]
            for j in range(k + 1, n + 1):
                row[j] = (row[j] * pivot - top[j] * factor) // prev
            row[k] = 0
        prev = pivot
    x = [None] * n
    for i in range(n - 1, -1, -1):
        diag = a[i][i]
        if diag == 0:
            return None
        acc = Fraction(a[i][n])
        row = a[i]
        for j in range(i + 1, n):
            acc -= row[j] * x[j]
        x[i] = acc / diag
    return x


_WINDOWS_35 = None


def _windows35():
    """Complete R(3,5,k) edge windows, loaded once from VALIDATION.json."""
    global _WINDOWS_35
    if _WINDOWS_35 is None:
        _WINDOWS_35 = r35_edge_windows(_DATA_DIR / VALIDATION_FILE)
    return _WINDOWS_35


def relaxed_rhs_intervals(order, edges):
    """Relaxed R2-R5 right-hand sides over surviving degree histograms.

    Mirrors the imported m=3 q-machinery's feasible set exactly: every
    integer degree histogram of the class whose per-neighborhood triangle
    window is nonempty (ceil/floor integrality on 3t).  R2 is
    class-determined and stays a point; R3/R4/R5 min/max over survivors,
    with R5 = t in [W_lo - q_hi, W_hi - q_lo] per the pinned convention.
    """
    order = _exact_int(order, "order")
    edges = _exact_int(edges, "edges")
    windows = _windows35()
    low, _high = degree_bounds(order)
    w_lo = w_hi = f_lo = f_hi = q_lo = q_hi = None
    for hist in degree_histograms(order, edges):
        wedge = tri3_lo = tri3_hi = fsum = 0
        for index, count in enumerate(hist):
            if not count:
                continue
            degree = low + index
            wedge += count * _c2(degree)
            tri3_lo += count * windows[degree][0]
            tri3_hi += count * windows[degree][1]
            fsum += count * _f_of(order, degree)
        tri_lo = -(-tri3_lo // 3)   # ceil; 3t >= tri3_lo, t integral
        tri_hi = tri3_hi // 3
        if tri_lo > tri_hi:
            continue
        w_lo = wedge if w_lo is None else min(w_lo, wedge)
        w_hi = wedge if w_hi is None else max(w_hi, wedge)
        f_lo = fsum if f_lo is None else min(f_lo, fsum)
        f_hi = fsum if f_hi is None else max(f_hi, fsum)
        q_lo = (wedge - tri_hi if q_lo is None
                else min(q_lo, wedge - tri_hi))
        q_hi = (wedge - tri_lo if q_hi is None
                else max(q_hi, wedge - tri_lo))
    if w_lo is None:
        raise InfeasibleEdgeClass(
            f"no degree histogram survives for ({order}, {edges})")
    exact_r2 = edges * _c2(order - 2)
    return {"R2": (exact_r2, exact_r2),
            "R3": (w_lo, w_hi),
            "R4": (f_lo, f_hi),
            "R5": (w_lo - q_hi, w_hi - q_lo)}


def enumerate_basic_solutions(order, edges):
    """Every feasible basis of the envelope LP, exactly and provably complete.

    Vertices of {s >= 0 : R1, R2 exact, R3-R5 two-sided} come from
    activating 9 of the 17 inequality rows.  Exactness-preserving cut: a
    basis activating both endpoints of one range row carries two identical
    coefficient rows and is singular; active > 3 forces that by pigeonhole
    over the 3 range rows.  Hence only zeroed_count in 6..9 with distinct
    range rows reach square solves; count and min/max equal the naive
    C(17, 9) = 24,310 enumeration exactly.  A zeroed pick coordinate is a
    legitimate vertex and is offered as Fraction(0).

    Returns (feasible_basis_count, {motif key: (Fraction min, Fraction max)}),
    where t is the tau-row over (order - 3), q is W(s) - t(s) with W from
    the R3 row, and diamond / k3k1 / i4 are raw coordinates.
    """
    order = _exact_int(order, "order")
    edges = _exact_int(edges, "edges")
    rhs = relaxed_rhs_intervals(order, edges)
    rows = relation_rows()
    r1 = rows["R1"]
    r2 = rows["R2"]
    eq1 = _c4(order)
    eq2 = rhs["R2"][0]
    # relaxed_rhs_intervals returns AGGREGATE intervals: R3 is the W range
    # and R5 the t range.  The LP rows are row-space: lift W through the
    # exact identity sum e^2 s = e*C(r-2,2) + 2W(r-4) + e(e-1), and lift t
    # through sum tau s = (r-3) t.  R4 is already row-space.
    w_lo, w_hi = rhs["R3"]
    r3_lo = eq2 + 2 * w_lo * (order - 4) + edges * (edges - 1)
    r3_hi = eq2 + 2 * w_hi * (order - 4) + edges * (edges - 1)
    t_lo, t_hi = rhs["R5"]
    ranges = ((rows["R3"], (r3_lo, r3_hi)),
              (rows["R4"], rhs["R4"]),
              (rows["R5"], ((order - 3) * t_lo, (order - 3) * t_hi)))
    count = 0
    best = dict.fromkeys(MOTIF4_KEYS)

    def offer(key, value):
        current = best[key]
        if current is None:
            best[key] = (value, value)
        else:
            lo, hi = current
            best[key] = (lo if lo <= value else value,
                         hi if hi >= value else value)

    for zeroed_count in range(6, 10):
        active = 9 - zeroed_count
        for zeroed in itertools.combinations(range(11), zeroed_count):
            keep = [v for v in range(11) if v not in zeroed]
            eq1_row = [r1[v] for v in keep] + [eq1]
            eq2_row = [r2[v] for v in keep] + [eq2]
            restricted = [([row[v] for v in keep], bounds)
                          for row, bounds in ranges]
            pick_positions = tuple(
                None if i in zeroed else keep.index(i)
                for i in _PICK_INDICES)
            for chosen_rows in itertools.combinations(range(3), active):
                for sides in itertools.product((0, 1), repeat=active):
                    mat = [eq1_row[:], eq2_row[:]]
                    for k, side in zip(chosen_rows, sides):
                        row_k, bounds = restricted[k]
                        mat.append(row_k + [bounds[side]])
                    sol = _solve_square_int(mat)
                    if sol is None or any(value < 0 for value in sol):
                        continue
                    feasible = True
                    vals = []
                    for row_k, bounds in restricted:
                        total = 0
                        for coef, x in zip(row_k, sol):
                            total += coef * x
                        if total < bounds[0] or total > bounds[1]:
                            feasible = False
                            break
                        vals.append(total)
                    if not feasible:
                        continue
                    count += 1
                    picked = [Fraction(0) if pos is None else sol[pos]
                              for pos in pick_positions]
                    offer("diamond", picked[0])
                    offer("k3k1", picked[1])
                    offer("i4", picked[2])
                    t_val = vals[2] / (order - 3)
                    offer("t", t_val)
                    w_val = (vals[0] - eq2 - edges * (edges - 1)) \
                        / (2 * (order - 4))
                    offer("q", w_val - t_val)
    return count, best


def verified_basic_solution_count(order, edges):
    """Exact feasible-basis count of the envelope LP for one class."""
    count, _best = enumerate_basic_solutions(order, edges)
    return count


_BASIC_CACHE = {}
_LP_FALLBACK = set()


def lp_fallback_classes():
    """Classes whose LP certification failed; consumed fail-closed."""
    return frozenset(_LP_FALLBACK)


def outer_motif_window(order, edges, key):
    """Certified outer window of one motif aggregate over an edge class.

    (floor(min), ceil(max)) over the envelope LP's verified basic solutions.
    ANY certification failure is fail-closed: the class is recorded in
    `lp_fallback_classes()` and the trivially sound (0, C(order, 4)) is
    returned.  `enumerate_basic_solutions` is resolved through the module
    namespace at call time (monkey-patchable by the audit suite).
    """
    if key not in MOTIF4_KEYS:
        raise ValueError(f"unknown motif aggregate {key!r}")
    cache_key = (order, edges)
    if cache_key not in _BASIC_CACHE:
        try:
            _BASIC_CACHE[cache_key] = enumerate_basic_solutions(order, edges)
        except Exception:
            _LP_FALLBACK.add(cache_key)
            _BASIC_CACHE[cache_key] = None
    entry = _BASIC_CACHE[cache_key]
    if entry is None:
        return (0, _c4(order))
    _count, best = entry
    window = best[key]
    if window is None:          # feasible histograms but empty polyhedron
        _LP_FALLBACK.add(cache_key)
        _BASIC_CACHE[cache_key] = None
        return (0, _c4(order))
    return (math.floor(window[0]), math.ceil(window[1]))


def catalog_motif_windows(records):
    """Exact per-class motif windows from streamed value-count histograms.

    `records`: {(order, edges): {key: {value: count}}}; fail-closed on any
    structural drift (missing motif keys, empty or negative histograms,
    nonpositive counts, or cross-key count drift).
    """
    if not records:
        raise ValueError("catalog motif windows need at least one class")
    windows = {}
    for cls, histograms in records.items():
        if (not isinstance(cls, tuple) or len(cls) != 2
                or type(cls[0]) is not int or type(cls[1]) is not int
                or cls[0] < 1 or not 0 <= cls[1] <= _c2(cls[0])):
            raise ValueError(f"malformed catalog class key {cls!r}")
        if not isinstance(histograms, dict) \
                or set(histograms) != set(MOTIF4_KEYS):
            raise ValueError(
                f"{cls}: histograms must cover exactly {MOTIF4_KEYS}")
        total = None
        merged = {}
        for key in MOTIF4_KEYS:
            histogram = histograms[key]
            if not histogram:
                raise ValueError(f"{cls}: empty {key} histogram")
            running = 0
            for value, count in histogram.items():
                if type(value) is not int:
                    raise ValueError(f"{cls}: non-integer {key} value {value!r}")
                if type(count) is not int or count <= 0:
                    raise ValueError(f"{cls}: nonpositive {key} count {count!r}")
                running += count
            if total is None:
                total = running
            elif running != total:
                raise ValueError(
                    f"{cls}: graph count drift across motif histograms")
            merged[key] = (min(histogram), max(histogram))
        windows[cls] = merged
    return windows


def motif_interval(order, edges, key, catalog):
    """Exact catalog window when present, else the certified outer window."""
    if catalog:
        window = catalog.get((order, edges))
        if window is not None:
            return window[key]
    return outer_motif_window(order, edges, key)


# =========================================================================
# The extended n=45 state cone
# =========================================================================


@dataclass(frozen=True, order=True)
class State:
    """One (degree, deficiency pair) cell of the n=45 cone, m=4-extended."""
    d: int
    a: int
    b: int
    deficiency: int
    excess_balance: int
    g_lo: int
    g_hi: int
    h_lo: int
    h_hi: int


_EDGE_BOUNDS = None
_M3_CATALOG_WINDOWS = None
_M3_INTERVAL_SOURCES = None
_M4_CATALOG_WINDOWS = None


def _edge_bounds_map():
    """Published R(4,5) edge extrema per order, loaded once."""
    global _EDGE_BOUNDS
    if _EDGE_BOUNDS is None:
        _EDGE_BOUNDS = load_edge_bounds(_DATA_DIR / TABLE_FILE)
    return _EDGE_BOUNDS


def _m3_catalog_windows():
    """Frozen m=3 catalog q-windows and interval sources, from the artifact.

    The m=3 artifact is read-only evidence; its own sweep-produced windows
    are exactly the ones its 3,215 g-records were built from, so reusing
    them reproduces the frozen domain by construction (test 7).
    """
    global _M3_CATALOG_WINDOWS, _M3_INTERVAL_SOURCES
    if _M3_CATALOG_WINDOWS is None:
        raw = _load_json(_M3_ARTIFACT_PATH)
        windows = {}
        for record in raw["catalog_windows"]:
            windows[(record["order"], record["edges"])] = (
                record["q_min"], record["q_max"])
        _M3_CATALOG_WINDOWS = windows
        _M3_INTERVAL_SOURCES = {
            (record["d"], record["a"], record["b"]):
                (record["x_interval_source"], record["y_interval_source"])
            for record in raw["states"]
        }
    return _M3_CATALOG_WINDOWS


def _m4_catalog_windows():
    """Exact catalog motif windows from the m=4 artifact, or {} if absent."""
    global _M4_CATALOG_WINDOWS
    if _M4_CATALOG_WINDOWS is None:
        if not _M4_ARTIFACT_PATH.exists():
            _M4_CATALOG_WINDOWS = {}
        else:
            raw = _load_json(_M4_ARTIFACT_PATH)
            _M4_CATALOG_WINDOWS = {
                (record["order"], record["edges"]): {
                    key: tuple(record["windows"][key]) for key in MOTIF4_KEYS
                }
                for record in raw["catalog_motif_windows"]
            }
    return _M4_CATALOG_WINDOWS


def _trivial_window(order):
    """The sound global motif window: no aggregate leaves [0, C(order, 4)]."""
    return (0, _c4(order))


def _resolve_motif(src_map, cls, key, order):
    """One motif window from a resolution map, trivial when the class is out."""
    windows = src_map.get(cls)
    if windows is not None and key in windows:
        return tuple(windows[key])
    return _trivial_window(order)


def _h_interval(d, a, b, motif_windows, motif_sources):
    """The m=4 h row interval of state (d, a, b), plain stratum classes.

    h(v) = 3*(47 - 2d)*t(x) + 4*diamond(x) - 6*k3k1(x) - 12*i4(y) with
    x = (d, E45(d) - a), y = (m, E45(m) - b), m = N_TARGET - 1 - d
    (DESIGN A: the y-side i4 resolves at the plain class, mirroring the
    frozen m=3 q(D_v) convention).
    """
    m = N_TARGET - 1 - d
    edge_bounds = _edge_bounds_map()
    x_edges = _edge_window(edge_bounds, d)[1] - a
    y_edges = _edge_window(edge_bounds, m)[1] - b
    x_cls = (d, x_edges)
    y_cls = (m, y_edges)
    parts = (
        affine_interval(3 * (N_TARGET + 2 - 2 * d),
                        _resolve_motif(motif_windows, x_cls, "t", d)),
        affine_interval(4, _resolve_motif(motif_windows, x_cls, "diamond", d)),
        affine_interval(-6, _resolve_motif(motif_windows, x_cls, "k3k1", d)),
        affine_interval(-12, _resolve_motif(motif_windows, y_cls, "i4", m)),
    )
    h_lo = sum(interval[0] for interval in parts)
    h_hi = sum(interval[1] for interval in parts)
    x_source = motif_sources.get(x_cls, TRIVIAL_SOURCE)
    y_source = motif_sources.get(y_cls, TRIVIAL_SOURCE)
    return (h_lo, h_hi), x_source, y_source


def _policy_motif_maps(motif_windows, motif_sources):
    """Default h-side policy: artifact windows, trivial source names."""
    if motif_windows is None:
        motif_windows = _m4_catalog_windows()
        motif_sources = {
            cls: MOTIF4_CATALOG_SOURCE for cls in motif_windows
        }
    elif motif_sources is None:
        motif_sources = {}
    return motif_windows, motif_sources


def m4_state_interval(d, a, b, motif_windows=None, motif_sources=None):
    """State (d, a, b) with the sound m=3 g interval and m=4 h interval."""
    base = m3_state_interval(
        d, a, b, _windows35(), _edge_bounds_map(), _m3_catalog_windows())
    motif_windows, motif_sources = _policy_motif_maps(
        motif_windows, motif_sources)
    (h_lo, h_hi), _x_src, _y_src = _h_interval(
        d, a, b, motif_windows, motif_sources)
    return State(d, a, b, base.deficiency, base.excess_balance,
                 base.g_lo, base.g_hi, h_lo, h_hi)


def build_states(motif_windows=None, motif_sources=None):
    """Every state of the n=45 cone, sorted; g side frozen-exact by import."""
    edge_bounds = _edge_bounds_map()
    m3_states = _m3_build_states(
        _windows35(), edge_bounds, _m3_catalog_windows())
    motif_windows, motif_sources = _policy_motif_maps(
        motif_windows, motif_sources)
    states = []
    for base in m3_states:
        (h_lo, h_hi), _x_src, _y_src = _h_interval(
            base.d, base.a, base.b, motif_windows, motif_sources)
        states.append(State(base.d, base.a, base.b, base.deficiency,
                            base.excess_balance, base.g_lo, base.g_hi,
                            h_lo, h_hi))
    return sorted(states)


def state_record(state, motif_windows=None, motif_sources=None):
    """One canonical v2 state record: nine fields plus four source strings."""
    _m3_catalog_windows()   # populates the interval-source map
    x_interval_source, y_interval_source = _M3_INTERVAL_SOURCES.get(
        (state.d, state.a, state.b), (OUTER_SOURCE, OUTER_SOURCE))
    motif_windows, motif_sources = _policy_motif_maps(
        motif_windows, motif_sources)
    m = N_TARGET - 1 - state.d
    edge_bounds = _edge_bounds_map()
    x_cls = (state.d, _edge_window(edge_bounds, state.d)[1] - state.a)
    y_cls = (m, _edge_window(edge_bounds, m)[1] - state.b)
    return {
        "d": state.d, "a": state.a, "b": state.b,
        "deficiency": state.deficiency,
        "excess_balance": state.excess_balance,
        "g_lo": state.g_lo, "g_hi": state.g_hi,
        "h_lo": state.h_lo, "h_hi": state.h_hi,
        "x_interval_source": x_interval_source,
        "y_interval_source": y_interval_source,
        "x_motif_source": motif_sources.get(x_cls, TRIVIAL_SOURCE),
        "y_motif_source": motif_sources.get(y_cls, TRIVIAL_SOURCE),
    }


# =========================================================================
# Execution pipeline (Task 2 Steps 4-5)
# =========================================================================

STREAM_BUDGET = 7200
ENVELOPE_BUDGET = 1800
EXECUTED_SOURCES = (
    "check_ramsey.py", "m3_deficiency_cone.py", "m4_deficiency_cone.py",
    "m4_subgraph_identities.py", "subgraph_identities.py",
)


def _stream_motif_counts(path, expected_order, expected_edges,
                         expected_count):
    """Exact 5-motif frequencies per (order, edges) for one graph6 catalog.

    One pass, no parsed graph retained; every line is order/edges-checked and
    the final line count must equal the validated census count.  Mirrors the
    frozen m=3 single-pass discipline with the m=4 kernel in place of q.
    """
    counts = {}
    graphs = 0
    try:
        stream = path.open("r", encoding="ascii")
    except OSError as exc:
        raise ValueError(f"{path.name}: unreadable catalog ({exc})") from None
    with stream:
        try:
            for number, line in enumerate(stream, start=1):
                text = line.strip()
                if not text:
                    continue
                graphs += 1
                try:
                    order, adj = parse_graph6_line(text)
                except ValueError as exc:
                    raise ValueError(
                        f"{path.name} line {number}: bad graph6 ({exc})"
                    ) from None
                if order != expected_order:
                    raise ValueError(
                        f"{path.name} line {number}: order {order}, expected "
                        f"{expected_order}"
                    )
                edges = edge_count(adj)
                if expected_edges is not None and edges != expected_edges:
                    raise ValueError(
                        f"{path.name} line {number}: {edges} edges, expected "
                        f"{expected_edges}"
                    )
                vector = motif4_vector(adj)
                frequency = counts.get((order, edges))
                if frequency is None:
                    frequency = counts[(order, edges)] = {
                        key: {} for key in MOTIF4_KEYS
                    }
                for key in MOTIF4_KEYS:
                    level = frequency[key]
                    value = vector[key]
                    level[value] = level.get(value, 0) + 1
        except (OSError, UnicodeError) as exc:
            raise ValueError(
                f"{path.name}: unreadable catalog ({exc})"
            ) from None
    if not graphs:
        raise ValueError(f"{path.name}: empty catalog")
    if graphs != expected_count:
        raise ValueError(
            f"{path.name}: {graphs} graphs, expected {expected_count}")
    return graphs, {
        cls: {key: {value: level[value] for value in sorted(level)}
              for key, level in sorted(freq.items())}
        for cls, freq in sorted(counts.items())
    }


def _replay_published_m4(path, expected_order, expected_graphs):
    """Zero m=4 specialized residual on every published graph and complement."""
    graphs = 0
    zeros = 0
    try:
        stream = path.open("r", encoding="ascii")
    except OSError as exc:
        raise ValueError(f"{path.name}: unreadable catalog ({exc})") from None
    with stream:
        try:
            for number, line in enumerate(stream, start=1):
                text = line.strip()
                if not text:
                    continue
                graphs += 1
                try:
                    order, adj = parse_graph6_line(text)
                except ValueError as exc:
                    raise ValueError(
                        f"{path.name} line {number}: bad graph6 ({exc})"
                    ) from None
                if order != expected_order:
                    raise ValueError(
                        f"{path.name} line {number}: order {order}, expected "
                        f"{expected_order}"
                    )
                for label, graph in (
                    ("graph", adj), ("complement", complement(adj)),
                ):
                    residual = m4_residual_specialized(graph)
                    if residual != 0:
                        raise ValueError(
                            f"{path.name} line {number}: {label} has residual "
                            f"{residual}, expected 0"
                        )
                    zeros += 1
        except (OSError, UnicodeError) as exc:
            raise ValueError(
                f"{path.name}: unreadable catalog ({exc})"
            ) from None
    if graphs != expected_graphs:
        raise ValueError(
            f"{path.name}: {graphs} graphs, expected {expected_graphs}")
    if zeros != 2 * expected_graphs:
        raise ValueError(
            f"{path.name}: {zeros} zero residuals, expected "
            f"{2 * expected_graphs}")
    return graphs, zeros


def _n49_replay(data_dir):
    """The MR97 Theorem 3.2 anchor on the two extremal (24, 132) graphs."""
    census_path = data_dir / R45_CENSUS_FILE
    extremal = []
    with census_path.open("r", encoding="ascii") as stream:
        for number, line in enumerate(stream, start=1):
            text = line.strip()
            if not text:
                continue
            order, adj = parse_graph6_line(text)
            if edge_count(adj) != 132:
                continue
            extremal.append(adj)
            if len(extremal) == 2:
                break
    if len(extremal) != 2:
        raise ValueError(
            f"{census_path.name}: expected two 132-edge extremal graphs, "
            f"found {len(extremal)}")
    observed = []
    for adj in extremal:
        if triangle_count(adj) != 176:
            raise ValueError("n=49 replay: triangles != 176")
        if diamond_count(adj) != 792:
            raise ValueError("n=49 replay: diamonds != 792")
        if induced_k3k1_count(adj) != 528:
            raise ValueError("n=49 replay: induced K3+K1 != 528")
        observed.append(independent_quad_count(adj))
    observed.sort()
    if observed != [138, 144]:
        raise ValueError(f"n=49 replay: observed i4 {observed} != [138, 144]")
    row = 3 * (49 + 2 - 2 * 24) * 176 + 4 * 792 - 6 * 528
    if row != 1584 or row != 12 * 132:
        raise ValueError(f"n=49 replay: row {row} != 1584 = 12 * 132")
    if min(observed) <= 132:
        raise ValueError("n=49 replay: extremal i4 not above the forced mean")
    return {"extremal_graphs": 2, "row_scaled": row, "forced_i4_mean": 132,
            "observed_i4": observed}


def _analyze(src_dir, echo):
    """The whole fail-closed m=4 analysis; returns the canonical document."""
    started = time.monotonic()
    data_dir = src_dir.parent / DATA_DIRNAME
    math_root = src_dir.parents[1]
    validation = data_dir / VALIDATION_FILE
    extreme_validation = data_dir / VALIDATION_EXTREME_FILE
    tables = data_dir / TABLE_FILE

    edge_bounds = load_edge_bounds(tables)
    windows35 = r35_edge_windows(validation)
    main_records = _select_main_records(validation)
    extreme_records = _select_extreme_records(extreme_validation)

    # Bind every consumed path before any expensive computation.
    inputs = []
    preflight = {}

    def add_input(path, graph_count, expected_sha=None):
        record = _input_record(path, math_root, graph_count, expected_sha)
        inputs.append(record)
        preflight[path] = record

    for name in R35_FILES:
        record = main_records[name]
        path = data_dir / name
        graphs = _count_graph6_lines(path)
        if graphs != record["graphs"]:
            raise ValueError(
                f"{name}: {graphs} graphs, expected {record['graphs']}")
        add_input(path, graphs, record["sha256"])
    for name in sorted(extreme_records):
        record = extreme_records[name]
        add_input(data_dir / R45_EXTREME_DIRNAME / name, record["graphs"])
    census = main_records[R45_CENSUS_FILE]
    census_path = data_dir / R45_CENSUS_FILE
    add_input(census_path, census["graphs"], census["sha256"])
    replay_path = data_dir / R55_REPLAY_FILE
    add_input(replay_path, main_records[R55_REPLAY_FILE]["graphs"],
              main_records[R55_REPLAY_FILE]["sha256"])
    for path in (validation, extreme_validation, tables):
        add_input(path, None)
    for name in EXECUTED_SOURCES:
        add_input(src_dir / name, None)
    # The frozen m=3 artifact anchors the g side of every state record.
    add_input(_M3_ARTIFACT_PATH, None)
    inputs = _sorted_inputs(inputs)

    replay_graphs, replay_zeros = _replay_published_m4(
        replay_path, R55_REPLAY_ORDER, R55_REPLAY_GRAPHS)
    echo(f"M4 REPLAY: {replay_zeros}/{2 * replay_graphs} "
         "specialized residual zero")

    n49 = _n49_replay(data_dir)
    echo("M4 N49 REPLAY: extremal=2 row=1584 forced_i4_mean=132 "
         f"observed={n49['observed_i4']}")

    extreme_names = sorted(extreme_records)
    catalog_records = [
        (data_dir / R45_EXTREME_DIRNAME / name,
         extreme_records[name]["order"],
         extreme_records[name]["edges"],
         extreme_records[name]["graphs"])
        for name in extreme_names
    ]
    catalog_records.append(
        (census_path, census["order"], None, census["graphs"]))

    histograms = {}
    graph_counts = {}
    owners = {}
    total_graphs = 0
    for path, order, edges, count in catalog_records:
        if time.monotonic() - started > STREAM_BUDGET:
            raise ValueError("catalog streaming exceeded the 7200 s budget")
        if path.as_posix() in graph_counts:
            raise ValueError(f"catalog {path.name} is listed twice")
        graphs, produced = _stream_motif_counts(path, order, edges, count)
        graph_counts[path.as_posix()] = graphs
        total_graphs += graphs
        for cls, frequency in produced.items():
            if cls in owners:
                raise ValueError(
                    f"class (order={cls[0]}, edges={cls[1]}) is produced by "
                    f"both {owners[cls].name} and {path.name}")
            owners[cls] = path
            histograms[cls] = frequency
    stream_seconds = format(time.monotonic() - started, ".1f")
    echo(f"M4 STREAM: {total_graphs} graphs in {stream_seconds} "
         f"(budget {STREAM_BUDGET})")

    catalog_windows = catalog_motif_windows(histograms)

    # Envelope self-audit: every observed value lies inside its LP window.
    envelope_start = time.monotonic()
    violations = 0
    for cls in sorted(catalog_windows):
        if time.monotonic() - envelope_start > ENVELOPE_BUDGET:
            raise ValueError("envelope LP exceeded the 1800 s budget")
        order, edges = cls
        for key in MOTIF4_KEYS:
            low, high = outer_motif_window(order, edges, key)
            observed = histograms[cls][key]
            violations += sum(
                1 for value in observed if not low <= value <= high)
    if violations:
        raise ValueError(
            f"{violations} exact catalog motif values fall outside the "
            "outer LP envelope")
    leaked = sorted(_LP_FALLBACK & set(catalog_windows))
    if leaked:
        raise ValueError(
            f"{len(leaked)} catalog classes silently fell back to trivia"
            f"l LP windows: {leaked}")
    audit_note = "FALLBACK" if _LP_FALLBACK else "ALL"
    echo(f"M4 ENVELOPES: catalog violations={violations} "
         f"outer-window audits={audit_note}")

    # Strata window resolution for the state cone.
    strata = set()
    edges_max = {order: bounds[1] for order, bounds in edge_bounds.items()}
    probes = build_states()     # trivial-policy pass; strata classes pinned
    for state in probes:
        m = N_TARGET - 1 - state.d
        strata.add((state.d, edges_max[state.d] - state.a))
        strata.add((m, edges_max[m] - state.b))
    motif_windows = {}
    motif_sources = {}
    outer_records = []
    for cls in sorted(strata):
        order, edges = cls
        if cls in catalog_windows:
            motif_windows[cls] = catalog_windows[cls]
            motif_sources[cls] = "catalog:" + _relative_path(
                owners[cls], math_root)
            method = "catalog"
            window_record = {key: list(catalog_windows[cls][key])
                             for key in MOTIF4_KEYS}
        else:
            if time.monotonic() - envelope_start > ENVELOPE_BUDGET:
                raise ValueError("envelope LP exceeded the 1800 s budget")
            resolved = {
                key: outer_motif_window(order, edges, key)
                for key in MOTIF4_KEYS
            }
            motif_windows[cls] = resolved
            if cls in _LP_FALLBACK:
                method, motif_sources[cls] = "trivial_fallback", TRIVIAL_SOURCE
            else:
                method, motif_sources[cls] = "envelope_lp", "envelope_lp"
            window_record = {key: list(resolved[key]) for key in MOTIF4_KEYS}
        outer_records.append({"order": order, "edges": edges,
                              "method": method, "windows": window_record})

    states = build_states(motif_windows, motif_sources)
    if not states:
        raise ValueError("the n=45 state cone is empty")
    for state in states:
        if state.h_lo > state.h_hi:
            raise ValueError(f"inverted h interval at {state}")
    echo(f"M4 STATES: {len(states)} states")

    for path, record in preflight.items():
        _verify_unchanged_input(path, record)

    def win_record(order):
        lo, hi = windows35[order]
        return {"order": order, "edge_min": lo, "edge_max": hi}

    document = {
        "schema_version": SCHEMA_VERSION,
        "campaign_id": CAMPAIGN_ID,
        "disposition": DISPOSITION,
        "inputs": inputs,
        "trust_roots": [
            {"statement": statement, "url": url}
            for statement, url in TRUST_ROOTS
        ],
        "replay": {"published_graphs": replay_graphs,
                   "complements": replay_graphs,
                   "residual_zero": replay_zeros},
        "n49_replay": n49,
        "r35_edge_windows": [win_record(order) for order in sorted(windows35)],
        "catalog_motif_histograms": [
            {"order": cls[0], "edges": cls[1],
             "source": _relative_path(owners[cls], math_root),
             "graph_count": sum(histograms[cls]["q"].values()),
             "histograms": {
                 key: [{"value": value, "count": freq[value]}
                       for value in sorted(freq)]
                 for key, freq in sorted(histograms[cls].items())
             }}
            for cls in sorted(histograms)
        ],
        "catalog_motif_windows": [
            {"order": cls[0], "edges": cls[1],
             "source": _relative_path(owners[cls], math_root),
             "windows": {key: list(catalog_windows[cls][key])
                         for key in MOTIF4_KEYS}}
            for cls in sorted(catalog_windows)
        ],
        "outer_motif_windows": outer_records,
        "states": [
            state_record(state, motif_windows, motif_sources)
            for state in states
        ],
        "searches": [],
    }
    return document


_TOP_LEVEL_KEYS = (
    "campaign_id", "catalog_motif_histograms", "catalog_motif_windows",
    "disposition", "inputs", "n49_replay", "outer_motif_windows",
    "r35_edge_windows", "replay", "schema_version", "states",
    "trust_roots", "searches",
)


def _validate_document(document):
    """Schema and count invariants of the v2 document, pre-write."""
    if sorted(document) != sorted(_TOP_LEVEL_KEYS):
        raise ValueError("v2 document top-level keys drifted")
    if document["schema_version"] != SCHEMA_VERSION:
        raise ValueError("schema version drifted")
    if document["campaign_id"] != CAMPAIGN_ID:
        raise ValueError("campaign id drifted")
    if document["disposition"] != DISPOSITION:
        raise ValueError("disposition drifted before Task 3")
    if document["searches"] != []:
        raise ValueError("Task 3 owns the searches section")
    if len(document["trust_roots"]) != 6:
        raise ValueError("exactly six trust roots are pinned")
    if not document["states"]:
        raise ValueError("the state cone must be nonempty")
    for record in document["states"]:
        if sorted(record) != sorted((
                "d", "a", "b", "deficiency", "excess_balance",
                "g_lo", "g_hi", "h_lo", "h_hi", "x_interval_source",
                "y_interval_source", "x_motif_source", "y_motif_source")):
            raise ValueError("state record keys drifted")
        for field in ("d", "a", "b", "deficiency", "excess_balance",
                      "g_lo", "g_hi", "h_lo", "h_hi"):
            if type(record[field]) is not int:
                raise ValueError(f"state field {field} is not an int")


def _write_document(target, document):
    """Serialize deterministically, parse back, validate v2, atomically replace."""
    text = json.dumps(document, indent=2, sort_keys=True, ensure_ascii=True)
    parent = target.parent
    if not parent.is_dir():
        raise ValueError(f"{target.name}: no output directory")
    handle, name = tempfile.mkstemp(
        dir=str(parent), prefix=target.name + ".", suffix=".tmp")
    temporary = Path(name)
    try:
        with os.fdopen(handle, "w", encoding="ascii") as stream:
            stream.write(text)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.chmod(name, 0o644)
        parsed = _load_json(temporary)
        if parsed != document:
            raise ValueError(
                "the serialized analysis does not parse back exactly")
        _validate_document(parsed)
        os.replace(name, str(target))
    except BaseException:
        try:
            temporary.unlink()
        except OSError:
            pass
        raise


def _echo(line):
    print(line, flush=True)


def main(argv=None):
    """Run the m=4 catalog analysis and write the v2 analysis JSON."""
    parser = argparse.ArgumentParser(
        description="Build the m=4 motif windows and the extended n=45 cone.")
    src_dir = _SRC_DIR
    default_output = src_dir.parent / DATA_DIRNAME / "higher_identity_m4.json"
    parser.add_argument(
        "--output", type=Path, default=default_output,
        help=f"analysis JSON target (default: {default_output})")
    args = parser.parse_args(argv)
    target = args.output
    if not target.is_absolute():
        target = (Path.cwd() / target).resolve()
    try:
        document = _analyze(src_dir, _echo)
        _validate_document(document)
        _write_document(target, document)
    except (ValueError, OSError) as exc:
        print(f"m4 analysis failed: {exc}", file=sys.stderr)
        return 1
    _echo(f"WROTE {_relative_path(target, src_dir.parents[1])}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
