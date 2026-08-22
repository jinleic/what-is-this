#!/usr/bin/env python3
"""Exact mixed-identity motif windows and the extended n=45 state cone.

This module is the v3 producer for the Engstrom K=4 mixed identity.  It
reuses the landed m4 relation system and state machinery for the g/h rows,
then adds the fifteen aggregate windows and the per-state F interval.
"""

from __future__ import annotations

import argparse
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

from check_ramsey import parse_graph6_line
from mixed_subgraph_identities import mixed_residual, mixed_vertex_row, motif_sweep
from subgraph_identities import complement, edge_count

import m4_deficiency_cone as _m4
from m3_deficiency_cone import (
    R35_FILES,
    R45_CENSUS_FILE,
    R45_EXTREME_DIRNAME,
    R55_REPLAY_FILE,
    R55_REPLAY_GRAPHS,
    R55_REPLAY_ORDER,
    DATA_DIRNAME,
    N_TARGET,
    OUTER_SOURCE,
    TABLE_FILE,
    VALIDATION_EXTREME_FILE,
    VALIDATION_FILE,
    _count_graph6_lines,
    _exact_int,
    _input_record,
    _load_json,
    _relative_path,
    _select_extreme_records,
    _select_main_records,
    _sorted_inputs,
    _verify_unchanged_input,
    degree_bounds,
    degree_histograms,
    load_edge_bounds,
    r35_edge_windows,
)


_SRC_DIR = Path(__file__).resolve().parent
_DATA_DIR = _SRC_DIR.parent / DATA_DIRNAME
_EXTREME_DIR = _DATA_DIR / R45_EXTREME_DIRNAME
_M3_ARTIFACT_PATH = _DATA_DIR / "higher_identity_m3.json"
_M4_ARTIFACT_PATH = _DATA_DIR / "higher_identity_m4.json"
_V3_ARTIFACT_PATH = _DATA_DIR / "engstrom_identity.json"

SCHEMA_VERSION = 3
CAMPAIGN_ID = "higher_order_identity_positive_deficiency_mixed"
DISPOSITION = "MIXED_ANALYZED"
EXPECTED_STATES = 3215
EXPECTED_CATALOG_GRAPHS = 8_500_211
TRIVIAL_SOURCE = "trivial:binomial_coefficient_bound"
DEFAULT_SOURCE = "mixed_default:m4_frozen_plus_trivial_new_keys"

# The order is frozen in the v3 schema and in every catalog histogram.
STREAM_KEYS = (
    "q", "t", "p3", "pc3", "i4", "k4", "diamond", "k3k1", "c4",
    "claw", "paw", "p4", "two_k2", "k2_2k1", "p3_k1",
)
MOTIF4_KEYS = tuple(_m4.MOTIF4_KEYS)
_M4_FROZEN5 = MOTIF4_KEYS

# Raw 4-vertex coordinates in relation_rows() order.
_RAW_KEYS = (
    "i4", "k2_2k1", "two_k2", "p3_k1", "p4", "claw", "k3k1",
    "c4", "paw", "diamond", "k4",
)
_RAW_INDEX = {key: index for index, key in enumerate(_RAW_KEYS)}

STREAM_BUDGET = 7200
ENVELOPE_BUDGET = 1800
EXECUTED_SOURCES = (
    "check_ramsey.py",
    "m3_deficiency_cone.py",
    "m4_deficiency_cone.py",
    "m4_subgraph_identities.py",
    "mixed_deficiency_cone.py",
    "mixed_subgraph_identities.py",
    "subgraph_identities.py",
)

TRUST_ROOTS = _m4.TRUST_ROOTS + (
    (
        "The Engstrom K=4 mixed identity itself is locally re-derived "
        "(module-docstring counting ledger, Task 1), exhaustively verified "
        "on all labeled graphs n<=6 (33,867 graphs), spot-verified at "
        "orders 7-8, replayed on all 656 n=42 graphs and complements, and "
        "replayed against the recorded n=49 type-combo values "
        "{0,144,288,432} with a unique zero corner on the frozen order-24 "
        "census. Primary-source anchors: Engstrom, arXiv:1002.4304, "
        "Theorem (Conjecture 1 in McKay and Radziszowski), and MR97 "
        "Conjecture 1.",
        "https://arxiv.org/abs/1002.4304",
    ),
    (
        "The local motif-window tables in data/engstrom_identity.json are "
        "derived only from hash-pinned frozen catalogs, the relation system "
        "R1-R5 (whose constants are unit-tested), and the frozen m=3 "
        "q-machinery; they are consumed downstream only through the "
        "artifact.",
        "https://arxiv.org/abs/2409.15709",
    ),
)


def _c2(value):
    return value * (value - 1) // 2


def _c3(value):
    return value * (value - 1) * (value - 2) // 6


def _c4(value):
    return value * (value - 1) * (value - 2) * (value - 3) // 24


def _sweep_vector(adj):
    """Return the canonical fifteen aggregate values for one graph."""
    sweep = motif_sweep(adj)
    if sweep.k4 != 0:
        raise ValueError("catalog graph contains a K4")
    return {
        "q": sweep.wedges - sweep.triangles,
        "t": sweep.triangles,
        "p3": sweep.induced_p3,
        "pc3": sweep.pc3,
        "i4": sweep.e4,
        "k4": sweep.k4,
        "diamond": sweep.diamond,
        "k3k1": sweep.k3k1,
        "c4": sweep.c4,
        "claw": sweep.claw,
        "paw": sweep.paw,
        "p4": sweep.p4,
        "two_k2": sweep.two_k2,
        "k2_2k1": sweep.k2_2k1,
        "p3_k1": sweep.p3_k1,
    }


def stream_catalog_motifs(
    path, expected_order=None, expected_edges=None, expected_count=None
):
    """Stream a catalog into fifteen histograms.

    The public one-argument form returns the histogram mapping from the plan;
    the producer/test form with explicit expectations returns ``(graphs,
    mapping)`` after enforcing those expectations.
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
                if expected_order is not None and order != expected_order:
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
                vector = _sweep_vector(adj)
                cls = (order, edges)
                frequency = counts.get(cls)
                if frequency is None:
                    frequency = counts[cls] = {key: {} for key in STREAM_KEYS}
                for key in STREAM_KEYS:
                    level = frequency[key]
                    value = vector[key]
                    level[value] = level.get(value, 0) + 1
        except (OSError, UnicodeError) as exc:
            raise ValueError(f"{path.name}: unreadable catalog ({exc})") from None
    if not graphs:
        raise ValueError(f"{path.name}: empty catalog")
    if expected_count is not None and graphs != expected_count:
        raise ValueError(
            f"{path.name}: {graphs} graphs, expected {expected_count}"
        )
    mapping = {
        cls: {
            key: {value: level[value] for value in sorted(level)}
            for key, level in sorted(freq.items())
        }
        for cls, freq in sorted(counts.items())
    }
    if expected_order is None and expected_edges is None and expected_count is None:
        return mapping
    return graphs, mapping


def catalog_motif_windows(records):
    """Convert value-count histograms to exact per-class min/max windows."""
    if not records:
        raise ValueError("catalog motif windows need at least one class")
    windows = {}
    for cls, histograms in records.items():
        if (
            not isinstance(cls, tuple)
            or len(cls) != 2
            or type(cls[0]) is not int
            or type(cls[1]) is not int
            or cls[0] < 1
            or not 0 <= cls[1] <= _c2(cls[0])
        ):
            raise ValueError(f"malformed catalog class key {cls!r}")
        if not isinstance(histograms, dict) or set(histograms) != set(STREAM_KEYS):
            raise ValueError(
                f"{cls}: histograms must cover exactly {STREAM_KEYS}"
            )
        total = None
        merged = {}
        for key in STREAM_KEYS:
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
                raise ValueError(f"{cls}: graph count drift across motif histograms")
            merged[key] = (min(histogram), max(histogram))
        windows[cls] = merged
    return windows


# ---------------------------------------------------------------------------
# Certified outer windows

_BASIC_CACHE = {}
_HISTOGRAM_CACHE = {}
_M4_OUTER_WINDOWS = None
_M4_CATALOG_WINDOWS = None
_M4_CATALOG_SOURCES = None
_V3_DEFAULT_CACHE = None
_LP_FALLBACK = set()


def _m4_outer_windows():
    """Read the frozen m4 outer table keyed by plain (order, edges)."""
    global _M4_OUTER_WINDOWS
    if _M4_OUTER_WINDOWS is None:
        if not _M4_ARTIFACT_PATH.exists():
            _M4_OUTER_WINDOWS = {}
        else:
            raw = _load_json(_M4_ARTIFACT_PATH)
            _M4_OUTER_WINDOWS = {}
            for record in raw.get("outer_motif_windows", []):
                cls = (record["order"], record["edges"])
                _M4_OUTER_WINDOWS[cls] = {
                    "method": record["method"],
                    "windows": {
                        key: tuple(record["windows"][key])
                        for key in MOTIF4_KEYS
                    },
                }
    return _M4_OUTER_WINDOWS


def _m4_catalog_windows():
    global _M4_CATALOG_WINDOWS, _M4_CATALOG_SOURCES
    if _M4_CATALOG_WINDOWS is None:
        if not _M4_ARTIFACT_PATH.exists():
            _M4_CATALOG_WINDOWS = {}
            _M4_CATALOG_SOURCES = {}
        else:
            raw = _load_json(_M4_ARTIFACT_PATH)
            _M4_CATALOG_WINDOWS = {}
            _M4_CATALOG_SOURCES = {}
            for record in raw.get("catalog_motif_windows", []):
                cls = (record["order"], record["edges"])
                _M4_CATALOG_WINDOWS[cls] = {
                    key: tuple(record["windows"][key])
                    for key in MOTIF4_KEYS
                }
                _M4_CATALOG_SOURCES[cls] = "catalog:" + record["source"]
    return _M4_CATALOG_WINDOWS


def _histogram_windows(order, edges):
    """Windows for q,t,p3,pc3 from surviving degree histograms."""
    cls = (order, edges)
    if cls in _HISTOGRAM_CACHE:
        return _HISTOGRAM_CACHE[cls]
    try:
        windows35 = _m4._windows35()
        low, _high = degree_bounds(order)
        result = {key: [None, None] for key in ("q", "t", "p3", "pc3")}
        for histogram in degree_histograms(order, edges):
            wedges = tri3_lo = tri3_hi = 0
            for index, count in enumerate(histogram):
                if not count:
                    continue
                degree = low + index
                wedges += count * _c2(degree)
                tri3_lo += count * windows35[degree][0]
                tri3_hi += count * windows35[degree][1]
            t_lo = -(-tri3_lo // 3)
            t_hi = tri3_hi // 3
            if t_lo > t_hi:
                continue
            values = {
                "q": (wedges - t_hi, wedges - t_lo),
                "t": (t_lo, t_hi),
                "p3": (wedges - 3 * t_hi, wedges - 3 * t_lo),
                "pc3": (
                    edges * (order - 2) - 2 * wedges + 3 * t_lo,
                    edges * (order - 2) - 2 * wedges + 3 * t_hi,
                ),
            }
            for key, (lo, hi) in values.items():
                current = result[key]
                current[0] = lo if current[0] is None else min(current[0], lo)
                current[1] = hi if current[1] is None else max(current[1], hi)
        if any(result[key][0] is None for key in result):
            raise ValueError(f"no degree histogram survives for {cls}")
        final = {key: tuple(value) for key, value in result.items()}
    except Exception:
        _LP_FALLBACK.add(cls)
        final = {key: (0, _c4(order)) for key in ("q", "t", "p3", "pc3")}
    _HISTOGRAM_CACHE[cls] = final
    return final


def _enumerate_basic_solutions(order, edges):
    """Enumerate the m4 envelope vertices while offering every raw coordinate."""
    order = _exact_int(order, "order")
    edges = _exact_int(edges, "edges")
    if order <= 4:
        raise ValueError("the relaxed envelope requires order > 4")
    rhs = _m4.relaxed_rhs_intervals(order, edges)
    rows = _m4.relation_rows()
    r1 = rows["R1"]
    r2 = rows["R2"]
    eq1 = _c4(order)
    eq2 = rhs["R2"][0]
    w_lo, w_hi = rhs["R3"]
    r3_lo = eq2 + 2 * w_lo * (order - 4) + edges * (edges - 1)
    r3_hi = eq2 + 2 * w_hi * (order - 4) + edges * (edges - 1)
    t_lo, t_hi = rhs["R5"]
    ranges = (
        (rows["R3"], (r3_lo, r3_hi)),
        (rows["R4"], rhs["R4"]),
        (rows["R5"], ((order - 3) * t_lo, (order - 3) * t_hi)),
    )
    count = 0
    best = dict.fromkeys(_RAW_KEYS + ("q", "t"))

    def offer(key, value):
        current = best[key]
        if current is None:
            best[key] = (value, value)
        else:
            lo, hi = current
            best[key] = (min(lo, value), max(hi, value))

    for zeroed_count in range(6, 10):
        active = 9 - zeroed_count
        for zeroed in itertools.combinations(range(11), zeroed_count):
            keep = [index for index in range(11) if index not in zeroed]
            eq1_row = [r1[index] for index in keep] + [eq1]
            eq2_row = [r2[index] for index in keep] + [eq2]
            restricted = [
                ([row[index] for index in keep], bounds)
                for row, bounds in ranges
            ]
            for chosen_rows in itertools.combinations(range(3), active):
                for sides in itertools.product((0, 1), repeat=active):
                    matrix = [eq1_row[:], eq2_row[:]]
                    for selected, side in zip(chosen_rows, sides):
                        row, bounds = restricted[selected]
                        matrix.append(row + [bounds[side]])
                    solution = _m4._solve_square_int(matrix)
                    if solution is None or any(value < 0 for value in solution):
                        continue
                    totals = []
                    feasible = True
                    for row, bounds in restricted:
                        total = sum(coef * value for coef, value in zip(row, solution))
                        if total < bounds[0] or total > bounds[1]:
                            feasible = False
                            break
                        totals.append(total)
                    if not feasible:
                        continue
                    count += 1
                    for key, index in _RAW_INDEX.items():
                        position = None if index in zeroed else keep.index(index)
                        offer(key, Fraction(0) if position is None else solution[position])
                    t_value = totals[2] / (order - 3)
                    w_value = (
                        totals[0] - eq2 - edges * (edges - 1)
                    ) / (2 * (order - 4))
                    offer("t", t_value)
                    offer("q", w_value - t_value)
    return count, best


def _lp_windows(order, edges):
    cls = (order, edges)
    if cls not in _BASIC_CACHE:
        try:
            _BASIC_CACHE[cls] = _enumerate_basic_solutions(order, edges)
        except Exception:
            _LP_FALLBACK.add(cls)
            _BASIC_CACHE[cls] = None
    entry = _BASIC_CACHE[cls]
    if entry is None:
        return {key: (0, _c4(order)) for key in _RAW_KEYS + ("q", "t")}
    _count, best = entry
    result = {}
    for key, value in best.items():
        if value is None:
            _LP_FALLBACK.add(cls)
            _BASIC_CACHE[cls] = None
            return {name: (0, _c4(order)) for name in _RAW_KEYS + ("q", "t")}
        result[key] = (math.floor(value[0]), math.ceil(value[1]))
    return result


def lp_fallback_classes():
    """Return classes whose exact outer certification fell back to trivial."""
    return frozenset(_LP_FALLBACK)


def outer_motif_window(order, edges, key):
    """Return a sound outer window for one of the fifteen aggregates."""
    if key not in STREAM_KEYS:
        raise ValueError(f"unknown motif aggregate {key!r}")
    cls = (order, edges)
    if key in MOTIF4_KEYS:
        frozen = _m4_outer_windows().get(cls)
        if frozen is not None:
            return frozen["windows"][key]
        return _lp_windows(order, edges)[key]
    if key in ("p3", "pc3"):
        return _histogram_windows(order, edges)[key]
    return _lp_windows(order, edges)[key]


def motif_interval(order, edges, key, catalog):
    """Use an exact catalog interval when supplied, otherwise an outer one."""
    if catalog:
        class_windows = catalog.get((order, edges))
        if class_windows is not None and key in class_windows:
            return tuple(class_windows[key])
    return outer_motif_window(order, edges, key)


# ---------------------------------------------------------------------------
# State cone and F interval

@dataclass(frozen=True, order=True)
class State:
    d: int
    a: int
    b: int
    deficiency: int
    excess_balance: int
    g_lo: int
    g_hi: int
    h_lo: int
    h_hi: int
    f_lo: int
    f_hi: int


def _edge_bounds_map():
    return _m4._edge_bounds_map()


def _m3_interval_sources():
    _m4._m3_catalog_windows()
    return _m4._M3_INTERVAL_SOURCES


def _trivial_window(order):
    return (0, _c4(order))


def _default_mixed_maps():
    """Build a cheap default map with frozen m4 h windows and sound F bounds."""
    global _V3_DEFAULT_CACHE
    if _V3_DEFAULT_CACHE is not None:
        return _V3_DEFAULT_CACHE
    windows = {}
    sources = {}
    m4_catalog = _m4_catalog_windows()
    for cls, record in _m4_outer_windows().items():
        base = m4_catalog.get(cls, record["windows"])
        windows[cls] = {
            key: tuple(base[key]) for key in MOTIF4_KEYS
        }
        for key in STREAM_KEYS:
            if key not in windows[cls]:
                windows[cls][key] = _trivial_window(cls[0])
        sources[cls] = DEFAULT_SOURCE
    _V3_DEFAULT_CACHE = windows, sources
    return _V3_DEFAULT_CACHE


def _policy_motif_maps(motif_windows, motif_sources):
    if motif_windows is None:
        motif_windows, default_sources = _default_mixed_maps()
        if motif_sources is None:
            motif_sources = default_sources
    elif motif_sources is None:
        motif_sources = {}
    return motif_windows, motif_sources


def _project_m4_windows(motif_windows):
    return {
        cls: {
            key: tuple(values[key])
            for key in MOTIF4_KEYS
            if key in values
        }
        for cls, values in motif_windows.items()
    }


def _f_interval(d, a, b, motif_windows):
    """Exact interval of the mixed F row over resolved motif windows."""
    n = N_TARGET
    m = n - 1 - d
    bounds = _edge_bounds_map()
    x_edges = _m4._edge_window(bounds, d)[1] - a
    z_edges = _m4._edge_window(bounds, m)[1] - b
    x_cls = (d, x_edges)
    z_cls = (m, z_edges)
    ex = x_edges
    ey = _c2(m) - z_edges
    exact = (
        n * (n - 3) * d
        - (n * n + 2 * n - 6) * d * d
        + 3 * n * d ** 3
        - 2 * d ** 4
        + 2 * (n * n + n - 8) * ex
        - 12 * ex * ex
        - 12 * (n - 1) * d * ex
        + 12 * d * d * ex
        + 4 * ey * ey
        - 2 * (n - 2) * d * ey
        + 4 * d * d * ey
    )
    terms = (
        (72, "c4", x_cls, d),
        (12 * (n - 2), "t", x_cls, d),
        (24, "claw", x_cls, d),
        (24, "p4", x_cls, d),
        (24, "paw", x_cls, d),
        (12 * (n + 2) - 24 * d, "p3", x_cls, d),
        (32, "diamond", x_cls, d),
        (2 * (n - 8) + 4 * d, "pc3", z_cls, m),
        (-12, "k3k1", z_cls, m),
        (-8, "two_k2", z_cls, m),
        (-8, "p3_k1", z_cls, m),
        (-24, "k2_2k1", z_cls, m),
    )
    intervals = [
        _m4.affine_interval(coef, motif_interval(order, edges, key, motif_windows))
        for coef, key, (order, edges), _order in terms
    ]
    return (
        exact + sum(interval[0] for interval in intervals),
        exact + sum(interval[1] for interval in intervals),
    )


def mixed_state_interval(d, a, b, motif_windows=None, motif_sources=None):
    motif_windows, motif_sources = _policy_motif_maps(
        motif_windows, motif_sources
    )
    m4_windows = _project_m4_windows(motif_windows)
    base = _m4.m4_state_interval(d, a, b, m4_windows, motif_sources)
    f_lo, f_hi = _f_interval(d, a, b, motif_windows)
    return State(
        d,
        a,
        b,
        base.deficiency,
        base.excess_balance,
        base.g_lo,
        base.g_hi,
        base.h_lo,
        base.h_hi,
        f_lo,
        f_hi,
    )


def build_states(motif_windows=None, motif_sources=None):
    """Build all 3,215 states, reusing m4 for every g/h endpoint."""
    motif_windows, motif_sources = _policy_motif_maps(
        motif_windows, motif_sources
    )
    m4_windows = _project_m4_windows(motif_windows)
    base_states = _m4.build_states(m4_windows, motif_sources)
    states = []
    for base in base_states:
        f_lo, f_hi = _f_interval(base.d, base.a, base.b, motif_windows)
        states.append(
            State(
                base.d,
                base.a,
                base.b,
                base.deficiency,
                base.excess_balance,
                base.g_lo,
                base.g_hi,
                base.h_lo,
                base.h_hi,
                f_lo,
                f_hi,
            )
        )
    return sorted(states)


def state_record(state, motif_windows=None, motif_sources=None):
    """Serialize one v3 state with interval and motif provenance."""
    motif_windows, motif_sources = _policy_motif_maps(
        motif_windows, motif_sources
    )
    sources = _m3_interval_sources()
    x_interval_source, y_interval_source = sources.get(
        (state.d, state.a, state.b), (OUTER_SOURCE, OUTER_SOURCE)
    )
    m = N_TARGET - 1 - state.d
    bounds = _edge_bounds_map()
    x_cls = (state.d, _m4._edge_window(bounds, state.d)[1] - state.a)
    y_cls = (m, _m4._edge_window(bounds, m)[1] - state.b)
    return {
        "d": state.d,
        "a": state.a,
        "b": state.b,
        "deficiency": state.deficiency,
        "excess_balance": state.excess_balance,
        "g_lo": state.g_lo,
        "g_hi": state.g_hi,
        "h_lo": state.h_lo,
        "h_hi": state.h_hi,
        "f_lo": state.f_lo,
        "f_hi": state.f_hi,
        "x_interval_source": x_interval_source,
        "y_interval_source": y_interval_source,
        "x_motif_source": motif_sources.get(x_cls, TRIVIAL_SOURCE),
        "y_motif_source": motif_sources.get(y_cls, TRIVIAL_SOURCE),
    }


# ---------------------------------------------------------------------------
# Replay and production analysis


def _replay_published_mixed(path, expected_order, expected_graphs):
    graphs = zeros = 0
    try:
        stream = path.open("r", encoding="ascii")
    except OSError as exc:
        raise ValueError(f"{path.name}: unreadable catalog ({exc})") from None
    with stream:
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
                    f"{path.name} line {number}: order {order}, expected {expected_order}"
                )
            for label, graph in (("graph", adj), ("complement", complement(adj))):
                residual = mixed_residual(graph)
                if residual != 0:
                    raise ValueError(
                        f"{path.name} line {number}: {label} residual {residual}"
                    )
                zeros += 1
    if graphs != expected_graphs:
        raise ValueError(f"{path.name}: {graphs} graphs, expected {expected_graphs}")
    if zeros != 2 * expected_graphs:
        raise ValueError(f"{path.name}: {zeros} zero residuals, expected {2 * expected_graphs}")
    return graphs, zeros


def _n49_replay(data_dir):
    census_path = data_dir / R45_CENSUS_FILE
    extremal = []
    with census_path.open("r", encoding="ascii") as stream:
        for line in stream:
            text = line.strip()
            if not text:
                continue
            order, adj = parse_graph6_line(text)
            if order == 24 and edge_count(adj) == 132:
                extremal.append(adj)
    if len(extremal) != 2:
        raise ValueError("n=49 replay: expected two 132-edge graphs")
    for adj in extremal:
        sweep = motif_sweep(adj)
        expected = {
            "triangles": 176,
            "diamond": 792,
            "paw": 1584,
            "claw": 792,
            "k3k1": 528,
        }
        for key, value in expected.items():
            observed = sweep.triangles if key == "triangles" else getattr(sweep, key)
            if observed != value:
                raise ValueError(f"n=49 replay: {key}={observed}, expected {value}")
    values = sorted(
        mixed_vertex_row(49, 24, motif_sweep(x), motif_sweep(z))
        for x in extremal for z in extremal
    )
    if values != [0, 144, 288, 432]:
        raise ValueError(f"n=49 replay: values {values} != [0, 144, 288, 432]")
    zero_corners = sum(value == 0 for value in values)
    if zero_corners != 1:
        raise ValueError(f"n=49 replay: zero corners {zero_corners}, expected 1")
    return {
        "extremal_graphs": 2,
        "combos": 4,
        "row_values": values,
        "zero_corners": zero_corners,
    }


def _strata_from_states(states, edge_bounds):
    edges_max = {order: bounds[1] for order, bounds in edge_bounds.items()}
    strata = set()
    for state in states:
        m = N_TARGET - 1 - state.d
        strata.add((state.d, edges_max[state.d] - state.a))
        strata.add((m, edges_max[m] - state.b))
    return strata


def _production_outer(
    cls, catalog_windows, histograms, owners, math_root, envelope_start=None
):
    """Resolve all fifteen windows for one stratum class."""
    if (
        envelope_start is not None
        and time.monotonic() - envelope_start > ENVELOPE_BUDGET
    ):
        raise ValueError("envelope certification exceeded the 1800 s budget")
    order, edges = cls
    if cls in catalog_windows:
        resolved = dict(catalog_windows[cls])
        return resolved, "catalog", "catalog:" + _relative_path(owners[cls], math_root)

    lp = _lp_windows(order, edges)
    histogram = _histogram_windows(order, edges)
    if (
        envelope_start is not None
        and time.monotonic() - envelope_start > ENVELOPE_BUDGET
    ):
        raise ValueError("envelope certification exceeded the 1800 s budget")
    frozen = _m4_outer_windows().get(cls)
    resolved = {}
    for key in STREAM_KEYS:
        if key in ("p3", "pc3"):
            resolved[key] = histogram[key]
        elif key in MOTIF4_KEYS:
            if frozen is not None:
                if key in lp and cls not in _LP_FALLBACK and lp[key] != frozen["windows"][key]:
                    raise ValueError(
                        f"m4 outer mismatch at {cls} {key}: {lp[key]} != "
                        f"{frozen['windows'][key]}"
                    )
                resolved[key] = frozen["windows"][key]
            else:
                resolved[key] = lp[key]
        else:
            resolved[key] = lp[key]
    if cls in _LP_FALLBACK:
        return (
            {key: _trivial_window(order) for key in STREAM_KEYS},
            "trivial_fallback",
            TRIVIAL_SOURCE,
        )
    return resolved, "envelope_lp", "envelope_lp"


def _analyze(src_dir, echo):
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
            raise ValueError(f"{name}: {graphs} graphs, expected {record['graphs']}")
        add_input(path, graphs, record["sha256"])
    for name in sorted(extreme_records):
        record = extreme_records[name]
        add_input(data_dir / R45_EXTREME_DIRNAME / name, record["graphs"])
    census = main_records[R45_CENSUS_FILE]
    census_path = data_dir / R45_CENSUS_FILE
    add_input(census_path, census["graphs"], census["sha256"])
    replay_path = data_dir / R55_REPLAY_FILE
    add_input(replay_path, main_records[R55_REPLAY_FILE]["graphs"], main_records[R55_REPLAY_FILE]["sha256"])
    for path in (validation, extreme_validation, tables, _M3_ARTIFACT_PATH, _M4_ARTIFACT_PATH):
        add_input(path, None)
    for name in EXECUTED_SOURCES:
        add_input(src_dir / name, None)
    inputs = _sorted_inputs(inputs)

    replay_graphs, replay_zeros = _replay_published_mixed(
        replay_path, R55_REPLAY_ORDER, R55_REPLAY_GRAPHS
    )
    echo(f"MIXED REPLAY: {replay_zeros}/{2 * replay_graphs} residual zero")
    n49 = _n49_replay(data_dir)
    echo(
        "MIXED N49 REPLAY: extremal=2 combos=4 "
        f"values={n49['row_values']} zero_corners={n49['zero_corners']}"
    )

    extreme_names = sorted(extreme_records)
    catalog_records = [
        (
            data_dir / R45_EXTREME_DIRNAME / name,
            extreme_records[name]["order"],
            extreme_records[name]["edges"],
            extreme_records[name]["graphs"],
        )
        for name in extreme_names
    ]
    catalog_records.append(
        (census_path, census["order"], None, census["graphs"])
    )

    histograms = {}
    owners = {}
    total_graphs = 0
    for path, order, edges, count in catalog_records:
        if time.monotonic() - started > STREAM_BUDGET:
            raise ValueError(f"catalog streaming exceeded the {STREAM_BUDGET} s budget")
        graphs, produced = stream_catalog_motifs(path, order, edges, count)
        if time.monotonic() - started > STREAM_BUDGET:
            raise ValueError(f"catalog streaming exceeded the {STREAM_BUDGET} s budget")
        total_graphs += graphs
        for cls, frequency in produced.items():
            if cls in owners:
                raise ValueError(f"class {cls} is produced by both catalogs")
            owners[cls] = path
            histograms[cls] = frequency
    expected_total = sum(count for _path, _order, _edges, count in catalog_records)
    if total_graphs != expected_total:
        raise ValueError(
            f"catalog total {total_graphs} graphs, expected {expected_total}"
        )
    if expected_total != EXPECTED_CATALOG_GRAPHS:
        raise ValueError(
            f"validated catalog total {expected_total}, expected "
            f"{EXPECTED_CATALOG_GRAPHS}"
        )
    echo(
        f"MIXED STREAM: {total_graphs} graphs in "
        f"{time.monotonic() - started:.1f}s (budget {STREAM_BUDGET})"
    )
    catalog_windows = catalog_motif_windows(histograms)

    # The five m4-shared aggregates must be byte-for-byte consistent.
    frozen_catalog = _m4_catalog_windows()
    for cls, windows in catalog_windows.items():
        if cls in frozen_catalog:
            for key in MOTIF4_KEYS:
                if windows[key] != frozen_catalog[cls][key]:
                    raise ValueError(f"shared catalog window mismatch at {cls} {key}")

    envelope_start = time.monotonic()
    violations = 0
    for cls, class_windows in sorted(catalog_windows.items()):
        if time.monotonic() - envelope_start > ENVELOPE_BUDGET:
            raise ValueError("envelope certification exceeded the 1800 s budget")
        order, edges = cls
        for key in STREAM_KEYS:
            low, high = outer_motif_window(order, edges, key)
            violations += sum(
                1
                for value in histograms[cls][key]
                if not low <= value <= high
            )
    if violations:
        raise ValueError(f"{violations} catalog values fall outside outer windows")
    leaked = sorted(_LP_FALLBACK & set(catalog_windows))
    if leaked:
        raise ValueError(f"catalog classes fell back to trivial certification: {leaked}")
    echo(
        f"MIXED ENVELOPES: catalog violations={violations} "
        f"outer-window audits={'FALLBACK' if _LP_FALLBACK else 'ALL'}"
    )

    # The state domain is independent of motif windows; the m4 artifact gives
    # the canonical projection used for the hard cross-check below.
    probes = build_states()
    strata = _strata_from_states(probes, edge_bounds)
    motif_windows = {}
    motif_sources = {}
    outer_records = []
    for cls in sorted(strata):
        resolved, method, source = _production_outer(
            cls, catalog_windows, histograms, owners, math_root, envelope_start
        )
        motif_windows[cls] = resolved
        motif_sources[cls] = source
        outer_records.append(
            {
                "order": cls[0],
                "edges": cls[1],
                "method": method,
                "windows": {key: list(resolved[key]) for key in STREAM_KEYS},
            }
        )

    if time.monotonic() - envelope_start > ENVELOPE_BUDGET:
        raise ValueError("envelope certification exceeded the 1800 s budget")
    states = build_states(motif_windows, motif_sources)
    if len(states) != EXPECTED_STATES:
        raise ValueError(f"mixed state count {len(states)}, expected {EXPECTED_STATES}")
    frozen_records = _load_json(_M4_ARTIFACT_PATH)["states"]
    if len(frozen_records) != EXPECTED_STATES:
        raise ValueError("frozen m4 state count drifted")
    frozen_states = {
        (record["d"], record["a"], record["b"]): record
        for record in frozen_records
    }
    if len(frozen_states) != EXPECTED_STATES:
        raise ValueError("frozen m4 state keys are not unique")
    for state in states:
        key = (state.d, state.a, state.b)
        record = frozen_states.get(key)
        if record is None:
            raise ValueError(f"mixed state missing from frozen m4 projection: {key}")
        for field in (
            "deficiency", "excess_balance", "g_lo", "g_hi", "h_lo", "h_hi"
        ):
            if getattr(state, field) != record[field]:
                raise ValueError(
                    f"mixed projection mismatch at {key} {field}: "
                    f"{getattr(state, field)} != {record[field]}"
                )
        rebuilt = state_record(state, motif_windows, motif_sources)
        for field in ("x_interval_source", "y_interval_source"):
            if rebuilt[field] != record[field]:
                raise ValueError(
                    f"mixed projection mismatch at {key} {field}: "
                    f"{rebuilt[field]!r} != {record[field]!r}"
                )
        if state.f_lo > state.f_hi:
            raise ValueError(f"inverted f interval at {state}")
    echo(f"MIXED STATES: {len(states)} states (asserted equal to frozen m4 projection)")

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
        "replay": {
            "published_graphs": replay_graphs,
            "complements": replay_graphs,
            "residual_zero": replay_zeros,
        },
        "n49_replay": n49,
        "r35_edge_windows": [win_record(order) for order in sorted(windows35)],
        "catalog_motif_histograms": [
            {
                "order": cls[0],
                "edges": cls[1],
                "source": _relative_path(owners[cls], math_root),
                "graph_count": sum(histograms[cls]["q"].values()),
                "histograms": {
                    key: [
                        {"value": value, "count": freq[value]}
                        for value in sorted(freq)
                    ]
                    for key, freq in sorted(histograms[cls].items())
                },
            }
            for cls in sorted(histograms)
        ],
        "catalog_motif_windows": [
            {
                "order": cls[0],
                "edges": cls[1],
                "source": _relative_path(owners[cls], math_root),
                "windows": {
                    key: list(catalog_windows[cls][key]) for key in STREAM_KEYS
                },
            }
            for cls in sorted(catalog_windows)
        ],
        "outer_motif_windows": outer_records,
        "states": [state_record(state, motif_windows, motif_sources) for state in states],
        "searches": [],
    }
    return document


_TOP_LEVEL_KEYS = (
    "campaign_id",
    "catalog_motif_histograms",
    "catalog_motif_windows",
    "disposition",
    "inputs",
    "n49_replay",
    "outer_motif_windows",
    "r35_edge_windows",
    "replay",
    "schema_version",
    "states",
    "trust_roots",
    "searches",
)


def _validate_document(document):
    if not isinstance(document, dict):
        raise ValueError("v3 document must be an object")
    if sorted(document) != sorted(_TOP_LEVEL_KEYS):
        raise ValueError("v3 document top-level keys drifted")
    if document["schema_version"] != SCHEMA_VERSION:
        raise ValueError("schema version drifted")
    if document["campaign_id"] != CAMPAIGN_ID:
        raise ValueError("campaign id drifted")
    if document["disposition"] != DISPOSITION:
        raise ValueError("disposition drifted before Task 3")
    if document["searches"] != []:
        raise ValueError("Task 3 owns the searches section")
    expected_roots = [
        {"statement": statement, "url": url}
        for statement, url in TRUST_ROOTS
    ]
    if document["trust_roots"] != expected_roots:
        raise ValueError("trust roots drifted")
    if document["replay"] != {
        "published_graphs": R55_REPLAY_GRAPHS,
        "complements": R55_REPLAY_GRAPHS,
        "residual_zero": 2 * R55_REPLAY_GRAPHS,
    }:
        raise ValueError("replay record drifted")
    if document["n49_replay"] != {
        "extremal_graphs": 2,
        "combos": 4,
        "row_values": [0, 144, 288, 432],
        "zero_corners": 1,
    }:
        raise ValueError("n49 replay record drifted")
    def valid_relative_path(value):
        if not isinstance(value, str) or not value or value.startswith("/"):
            return False
        parts = Path(value).parts
        return ".." not in parts and "." not in parts

    inputs = document["inputs"]
    if not isinstance(inputs, list):
        raise ValueError("inputs must be a list")
    input_paths = []
    for record in inputs:
        if set(record) != {"relative_path", "sha256", "bytes", "graph_count"}:
            raise ValueError("input record keys drifted")
        digest = record["sha256"]
        if (
            not isinstance(digest, str)
            or len(digest) != 64
            or any(char not in "0123456789abcdef" for char in digest)
        ):
            raise ValueError("input digest is not a lowercase SHA-256")
        if (
            not valid_relative_path(record["relative_path"])
            or type(record["bytes"]) is not int
            or record["bytes"] <= 0
        ):
            raise ValueError("malformed input record")
        count = record["graph_count"]
        if count is not None and (type(count) is not int or count <= 0):
            raise ValueError("malformed input graph count")
        input_paths.append(record["relative_path"])
    if input_paths != sorted(input_paths) or len(input_paths) != len(set(input_paths)):
        raise ValueError("inputs are not sorted and unique")
    windows = document["r35_edge_windows"]
    if not isinstance(windows, list):
        raise ValueError("r35 edge windows must be a list")
    previous = None
    for record in windows:
        if set(record) != {"order", "edge_min", "edge_max"}:
            raise ValueError("r35 edge-window keys drifted")
        order = record["order"]
        if (
            type(order) is not int
            or type(record["edge_min"]) is not int
            or type(record["edge_max"]) is not int
            or record["edge_min"] > record["edge_max"]
            or (previous is not None and order <= previous)
        ):
            raise ValueError("r35 edge windows are malformed or unsorted")
        previous = order

    def validate_class_records(records, expected_keys, label):
        if not isinstance(records, list):
            raise ValueError(f"{label} must be a list")
        seen = set()
        for record in records:
            if label == "outer_motif_windows":
                expected_record_keys = {"order", "edges", "method", "windows"}
            elif label == "catalog_motif_windows":
                expected_record_keys = {"order", "edges", "source", "windows"}
            else:
                expected_record_keys = {
                    "order", "edges", "source", "graph_count", "histograms"
                }
            if set(record) != expected_record_keys:
                raise ValueError(f"{label} record keys drifted")
            cls = (record["order"], record["edges"])
            if (
                type(cls[0]) is not int
                or type(cls[1]) is not int
                or cls in seen
                or (seen and cls <= max(seen))
            ):
                raise ValueError(f"{label} classes are malformed or unsorted")
            seen.add(cls)
            if label != "outer_motif_windows" and not valid_relative_path(
                record["source"]
            ):
                raise ValueError(f"{label} source is not relative")
            if label == "catalog_motif_histograms":
                if type(record["graph_count"]) is not int or record["graph_count"] <= 0:
                    raise ValueError("catalog graph count is malformed")
                histograms = record["histograms"]
                if set(histograms) != set(STREAM_KEYS):
                    raise ValueError("catalog histogram keys drifted")
                for entries in histograms.values():
                    if not isinstance(entries, list) or not entries:
                        raise ValueError("empty catalog histogram")
                    total = 0
                    values = []
                    for entry in entries:
                        if (
                            set(entry) != {"value", "count"}
                            or type(entry["value"]) is not int
                            or type(entry["count"]) is not int
                            or entry["count"] <= 0
                        ):
                            raise ValueError("malformed catalog histogram entry")
                        values.append(entry["value"])
                        total += entry["count"]
                    if values != sorted(values) or total != record["graph_count"]:
                        raise ValueError("catalog histogram ordering/count drifted")
            else:
                if set(record["windows"]) != set(STREAM_KEYS):
                    raise ValueError(f"{label} window keys drifted")
                for lo, hi in record["windows"].values():
                    if (
                        type(lo) is not int
                        or type(hi) is not int
                        or lo > hi
                    ):
                        raise ValueError(f"malformed {label} interval")
                if label == "outer_motif_windows" and record["method"] not in {
                    "catalog", "envelope_lp", "trivial_fallback"
                }:
                    raise ValueError("unknown outer-window method")
        return seen

    catalog_records = validate_class_records(
        document["catalog_motif_histograms"], STREAM_KEYS,
        "catalog_motif_histograms",
    )
    catalog_windows = validate_class_records(
        document["catalog_motif_windows"], STREAM_KEYS,
        "catalog_motif_windows",
    )
    outer_records = validate_class_records(
        document["outer_motif_windows"], STREAM_KEYS,
        "outer_motif_windows",
    )
    if catalog_records != catalog_windows:
        raise ValueError("catalog histogram/window class sets drifted")
    histogram_windows = {}
    for record in document["catalog_motif_histograms"]:
        cls = (record["order"], record["edges"])
        histogram_windows[cls] = {
            key: (entries[0]["value"], entries[-1]["value"])
            for key, entries in record["histograms"].items()
        }
    claimed_windows = {
        (record["order"], record["edges"]): {
            key: tuple(record["windows"][key])
            for key in STREAM_KEYS
        }
        for record in document["catalog_motif_windows"]
    }
    if histogram_windows != claimed_windows:
        raise ValueError("catalog windows do not match their histograms")
    if outer_records != set(_m4_outer_windows()):
        raise ValueError("outer-window class set drifted from frozen m4 strata")
    states = document["states"]
    if len(states) != EXPECTED_STATES:
        raise ValueError("state count drifted")
    state_keys = []
    expected_state_fields = {
        "d", "a", "b", "deficiency", "excess_balance", "g_lo", "g_hi",
        "h_lo", "h_hi", "f_lo", "f_hi", "x_interval_source",
        "y_interval_source", "x_motif_source", "y_motif_source",
    }
    for record in states:
        if set(record) != expected_state_fields:
            raise ValueError("state record keys drifted")
        for field in (
            "d", "a", "b", "deficiency", "excess_balance", "g_lo", "g_hi",
            "h_lo", "h_hi", "f_lo", "f_hi",
        ):
            if type(record[field]) is not int:
                raise ValueError(f"state field {field} is not an int")
        for field in (
            "x_interval_source", "y_interval_source",
            "x_motif_source", "y_motif_source",
        ):
            if not isinstance(record[field], str) or not record[field]:
                raise ValueError(f"state source {field} is malformed")
        if record["g_lo"] > record["g_hi"]:
            raise ValueError("inverted g interval")
        if record["h_lo"] > record["h_hi"]:
            raise ValueError("inverted h interval")
        if record["f_lo"] > record["f_hi"]:
            raise ValueError("inverted f interval")
        state_keys.append((record["d"], record["a"], record["b"]))
    if state_keys != sorted(state_keys) or len(state_keys) != len(set(state_keys)):
        raise ValueError("states are not sorted and unique")


def _write_document(target, document):
    text = json.dumps(document, indent=2, sort_keys=True, ensure_ascii=True)
    parent = target.parent
    if not parent.is_dir():
        raise ValueError(f"{target.name}: no output directory")
    handle, name = tempfile.mkstemp(
        dir=str(parent), prefix=target.name + ".", suffix=".tmp"
    )
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
            raise ValueError("serialized analysis does not parse back exactly")
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
    parser = argparse.ArgumentParser(
        description="Build the mixed motif windows and extended n=45 cone."
    )
    default_output = _DATA_DIR / "engstrom_identity.json"
    parser.add_argument("--output", type=Path, default=default_output)
    args = parser.parse_args(argv)
    target = args.output
    if not target.is_absolute():
        target = (Path.cwd() / target).resolve()
    try:
        document = _analyze(_SRC_DIR, _echo)
        _validate_document(document)
        _write_document(target, document)
    except (OSError, ValueError) as exc:
        print(f"mixed analysis failed: {exc}", file=sys.stderr)
        return 1
    _echo(f"WROTE {_relative_path(target, _SRC_DIR.parents[1])}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
