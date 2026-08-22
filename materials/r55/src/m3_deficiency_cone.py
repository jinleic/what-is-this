#!/usr/bin/env python3
"""Sound motif envelopes and the finite n=45 deficiency state cone (m=3 core).

Setting
-------
Suppose G is an R(5,5,45) graph. For a vertex v of degree d write
X = G_v^+ for the neighborhood and K for the complement of the
non-neighborhood G_v^-. Then X is an R(4,5,d) graph, K is an
R(4,5,m) graph with m = 44 - d, and the R(4,5) order bound R(4,5) = 25
forces d <= 24 and m <= 24, i.e. d in {20, ..., 24}.

The deficiency coordinates of the vertex are the two edge shortfalls against
the published R(4,5) edge maxima E_45:

    e(X) = E_45(d) - a,        e(K) = E_45(m) - b,        a, b >= 0.

Both shortfalls are bounded by the published edge windows, so
0 <= a <= E_45(d) - e_45(d) and 0 <= b <= E_45(m) - e_45(m). The pair
(d, a, b) is a *state*; the whole n=45 cone is a finite set of states.

The m=3 contribution of a state
-------------------------------
`subgraph_identities` proves, for every graph and every vertex,

    3 t(G_v^-) = (n + 3 - 3 d) e(X) + 3 p_3(X) + 6 t(X).

With q(H) = sum_x binom(d_x, 2) - t(H) and p_3 = q - 2t one has
3 p_3(X) + 6 t(X) = 3 q(X), and t(G_v^-) is the number of independent
triples of K, so the four-vertex identity

    i_3(H) = binom(r, 3) - e (r - 2) + q(H)          (H of order r, e edges)

turns the vertex identity into a statement about two R(4,5) graphs only:

    g_v := 3 [ binom(m, 3) - e(K) (m - 2) + q(K) ]
           - (48 - 3 d) e(X) - 3 q(X)  =  0.

Every term is an exact integer. A state pins the orders and the edge counts
but not q, so g_v is bracketed by substituting an interval for each q. The
bracket is what this module computes: `g_lo <= g_v <= g_hi`.

Sound q envelopes
-----------------
For an R(4,5,r,e) graph H the neighborhood of each vertex is K3-free and
I5-free, hence an R(3,5) graph, and H minus a closed neighborhood is K4-free
and I4-free, hence an R(4,4) graph. With R(3,5) = 14 and R(4,4) = 18 this
gives the degree window max(0, r - 18) <= d_x <= min(13, r - 1), and the
complete R(3,5,k) edge tables give

    sum_x e_35^min(d_x) <= 3 t(H) <= sum_x e_35^max(d_x),

because summing e(N_H(x)) over x counts each triangle of H three times.
`q_outer_interval` therefore enumerates *every integer* degree histogram that
satisfies the degree window and the handshake equation - graphicality is
deliberately not imposed, which only enlarges the family and keeps the
envelope an outer relaxation - and takes the extreme values of
wedges - t over that family. Integrality is used once more: 3t lies in an
integer interval, so t lies in [ceil(lo/3), floor(hi/3)].

Where a catalog is explicitly claimed complete by the Angeltveit-McKay data,
the exact observed min/max of q replace the outer interval for that
(order, edges) class. Missing classes keep the outer interval, so an unlisted
positive-deficiency class is never assumed empty. The published completeness
of the order-24 census and of the edge-extremal catalogs is a named trust
root: local hashing and validation do not re-prove it.

Link to the m=2 excess identity
-------------------------------
`structural_constraints` uses the exact per-vertex excess
e(G_v^-) - e(X) - d (n - 2 d) / 2, which sums to zero over all vertices.
Substituting the two deficiency coordinates and scaling by two gives

    excess_balance(v) = 2 (a + b) - budget2(d),
    budget2(d) = d (45 - 2 d) - 2 [ binom(m, 2) - E_45(m) - E_45(d) ],

so sum_v excess_balance(v) = 0. Together with sum_v g_v = 0 these are the
only global relations the affine certificate of Task 3 may use.

`data/structural_tables.json` is the single owner of the R(4,5) edge extrema;
this module loads them and derives budget2 instead of copying constants.

Calibration CLI
---------------
`main` selects only fully validated catalog records, hashes every consumed
file, replays the m=3 identity on the 328 published R(5,5,42) graphs and
their complements, streams the validated R(4,5) catalogs into exact q
histograms, checks every observed q against `q_outer_interval`, and writes
the canonical analysis document to `--output`. That document carries no
timestamp and no absolute path - every path is a canonical `math/`-relative
string - so it is reproducible byte for byte from the same inputs. It is
written to a temporary sibling, parsed back, and schema-checked before an
atomic replace, so a failed run never damages a good artifact.
"""

import argparse
import hashlib
import json
import os
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_ramsey import parse_graph6_line, popcount  # noqa: E402
from subgraph_identities import (  # noqa: E402
    complement, m3_residual_scaled, triangle_count,
)

N_TARGET = 45        # order of the hypothetical R(5,5,45) graph
R44 = 18             # R(4,4): bounds a non-neighborhood inside an R(4,5) graph
R35_MAX_ORDER = 13   # R(3,5) - 1: maximum degree of an R(4,5) graph
R45_MAX_ORDER = 24   # R(4,5) - 1: maximum order of an R(4,5) graph
# Both X and K are R(4,5) graphs, so d <= 24 and N_TARGET - 1 - d <= 24.
STATE_DEGREES = tuple(
    range(N_TARGET - 1 - R45_MAX_ORDER, R45_MAX_ORDER + 1)
)


class InfeasibleEdgeClass(ValueError):
    """No integer degree histogram survives for an (order, edges) class."""


@dataclass(frozen=True, order=True)
class State:
    """One (degree, deficiency pair) cell of the n=45 cone."""

    d: int
    a: int
    b: int
    deficiency: int
    excess_balance: int
    g_lo: int
    g_hi: int


def _exact_int(value, label):
    """Reject anything that is not already an exact Python int."""
    if type(value) is not int:
        raise ValueError(f"{label} must be an exact integer, got {value!r}")
    return value


def _load_json(path):
    try:
        raw = json.loads(path.read_text(encoding="ascii"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"{path}: unreadable JSON ({exc})") from None
    if not isinstance(raw, dict):
        raise ValueError(f"{path}: top level JSON object expected")
    return raw


def _window(record, low_key, high_key, label):
    low = _exact_int(record.get(low_key), f"{label} {low_key}")
    high = _exact_int(record.get(high_key), f"{label} {high_key}")
    if low > high:
        raise ValueError(f"{label}: empty window ({low} > {high})")
    return low, high


def load_edge_bounds(path):
    """Published R(4,5) edge extrema per order, from the owning table."""
    raw = _load_json(path)
    table = raw.get("extremal_edge_tables_R45")
    if not isinstance(table, dict):
        raise ValueError(f"{path}: no extremal_edge_tables_R45 object")
    bounds = {}
    for order, values in table.items():
        if not isinstance(values, dict):
            raise ValueError(f"{path}: order {order!r} is not an object")
        bounds[int(order)] = _window(values, "min", "max", f"order {order}")
    missing = [
        order for order in _cone_orders() if order not in bounds
    ]
    if missing:
        raise ValueError(f"{path}: no edge extrema for orders {missing}")
    return bounds


def _cone_orders():
    """Every R(4,5) order the n=45 cone needs: the degrees and their duals."""
    return sorted(
        set(STATE_DEGREES) | {N_TARGET - 1 - d for d in STATE_DEGREES}
    )


def r35_edge_windows(validation_path):
    """Complete R(3,5,k) edge windows, order zero synthesized as (0, 0)."""
    raw = _load_json(validation_path)
    results = raw.get("results")
    if not isinstance(results, list):
        raise ValueError(f"{validation_path}: no results list")
    windows = {0: (0, 0)}  # the empty graph: no vertices, no edges
    for record in results:
        if not isinstance(record, dict):
            raise ValueError(f"{validation_path}: malformed result record")
        if (record.get("s"), record.get("t")) != (3, 5):
            continue
        if not (record.get("ok") and record.get("count_matches")):
            continue  # unvalidated records are dropped, then reported missing
        orders = record.get("n_values")
        if not isinstance(orders, dict) or len(orders) != 1:
            raise ValueError(
                f"{validation_path}: record {record.get('file')!r} does not "
                "cover exactly one order"
            )
        order = int(next(iter(orders)))
        window = _window(record, "edge_min", "edge_max", f"R(3,5,{order})")
        if windows.setdefault(order, window) != window:
            raise ValueError(
                f"{validation_path}: conflicting windows for R(3,5,{order})"
            )
    missing = [
        order for order in range(R35_MAX_ORDER + 1) if order not in windows
    ]
    if missing:
        raise ValueError(
            f"{validation_path}: incomplete R(3,5) tables, missing {missing}"
        )
    return windows


def degree_bounds(order):
    """Ramsey degree window of an R(4,5,order) graph."""
    return max(0, order - R44), min(R35_MAX_ORDER, order - 1)


def degree_histograms(order, edges):
    """Every integer degree histogram of an R(4,5,order,edges) graph.

    Counts are indexed from `degree_bounds(order)[0]`. Graphicality is not
    imposed; only the degree window and the handshake equation are, so the
    family is an outer relaxation of the realizable histograms. Yielded in
    lexicographic order.
    """
    order = _exact_int(order, "order")
    edges = _exact_int(edges, "edges")
    if order < 0 or edges < 0:
        raise ValueError(f"negative class (order={order}, edges={edges})")
    low, high = degree_bounds(order)
    return _walk_histograms(low, high, order, 2 * edges)


def _walk_histograms(low, high, vertices, degree_sum):
    width = high - low + 1
    if width <= 0:
        return
    counts = [0] * width

    def walk(index, vertices_left, sum_left):
        degree = low + index
        if index == width - 1:
            if degree * vertices_left == sum_left:
                counts[index] = vertices_left
                yield tuple(counts)
                counts[index] = 0
            return
        for count in range(vertices_left + 1):
            sum_rest = sum_left - count * degree
            if sum_rest < 0:
                break
            vertices_rest = vertices_left - count
            if sum_rest > vertices_rest * high:
                break  # a larger count only widens this deficit
            if sum_rest < vertices_rest * (degree + 1):
                continue  # too many vertices left for the remaining degrees
            counts[index] = count
            yield from walk(index + 1, vertices_rest, sum_rest)
            counts[index] = 0

    yield from walk(0, vertices, degree_sum)


def q_outer_interval(order, edges, windows):
    """Outer interval for q = sum_x binom(d_x, 2) - t over an edge class."""
    low, high = degree_bounds(order)
    degrees = range(low, high + 1)
    missing = [degree for degree in degrees if degree not in windows]
    if missing:
        raise ValueError(
            f"no R(3,5) edge window for degrees {missing} "
            f"(order={order}, edges={edges})"
        )
    wedges = [degree * (degree - 1) // 2 for degree in degrees]
    lows = [_exact_int(windows[degree][0], f"e35_min({degree})")
            for degree in degrees]
    highs = [_exact_int(windows[degree][1], f"e35_max({degree})")
             for degree in degrees]
    q_lo = None
    q_hi = None
    histograms = 0
    for hist in degree_histograms(order, edges):
        histograms += 1
        wedge_sum = 0
        tri3_lo = 0
        tri3_hi = 0
        for count, wedge, low_edges, high_edges in zip(
            hist, wedges, lows, highs
        ):
            if count:
                wedge_sum += count * wedge
                tri3_lo += count * low_edges
                tri3_hi += count * high_edges
        triangle_lo = (tri3_lo + 2) // 3  # 3t >= tri3_lo, t integral
        triangle_hi = tri3_hi // 3
        if triangle_lo > triangle_hi:
            continue
        local_lo = wedge_sum - triangle_hi
        local_hi = wedge_sum - triangle_lo
        q_lo = local_lo if q_lo is None else min(q_lo, local_lo)
        q_hi = local_hi if q_hi is None else max(q_hi, local_hi)
    if q_lo is None:
        if histograms:
            raise InfeasibleEdgeClass(
                f"no histogram meets the R(3,5) triangle window for "
                f"order={order}, edges={edges}"
            )
        raise InfeasibleEdgeClass(
            f"no degree histogram for order={order}, edges={edges}"
        )
    return q_lo, q_hi


def q_interval(order, edges, windows, catalog_windows):
    """Exact catalog window when the class is complete, else the relaxation."""
    key = (order, edges)
    if key in catalog_windows:
        low, high = catalog_windows[key]
        low = _exact_int(low, f"catalog q_min{key}")
        high = _exact_int(high, f"catalog q_max{key}")
        if low > high:
            raise ValueError(f"catalog window {key} is empty ({low} > {high})")
        return low, high
    return q_outer_interval(order, edges, windows)


def budget2(d, edge_bounds):
    """Twice the largest m=2 excess deficit a degree-d vertex can carry."""
    m = N_TARGET - 1 - d
    maxima = _edge_maxima(edge_bounds, (d, m))
    h2_min = (
        2 * (m * (m - 1) // 2 - maxima[m] - maxima[d])
        - d * (N_TARGET - 2 * d)
    )
    return -h2_min


def _edge_maxima(edge_bounds, orders):
    return {order: _edge_window(edge_bounds, order)[1] for order in orders}


def _edge_window(edge_bounds, order):
    try:
        low, high = edge_bounds[order]
    except KeyError:
        raise ValueError(f"no R(4,5) edge extrema for order {order}") from None
    low = _exact_int(low, f"e45({order})")
    high = _exact_int(high, f"E45({order})")
    if low > high:
        raise ValueError(f"empty edge window for order {order}")
    return low, high


def _state_from_q(d, a, b, ex, ey, x_interval, y_interval, budget):
    """Assemble one state from the two q intervals; the sole g formula."""
    m = N_TARGET - 1 - d
    qx_lo, qx_hi = x_interval
    qy_lo, qy_hi = y_interval
    dual_base = 3 * (m * (m - 1) * (m - 2) // 6 - ey * (m - 2))
    rhs_base = (N_TARGET + 3 - 3 * d) * ex
    return State(
        d=d,
        a=a,
        b=b,
        deficiency=a + b,
        excess_balance=2 * (a + b) - budget,
        g_lo=dual_base + 3 * qy_lo - rhs_base - 3 * qx_hi,
        g_hi=dual_base + 3 * qy_hi - rhs_base - 3 * qx_lo,
    )


def m3_state_interval(d, a, b, windows, edge_bounds, catalog_windows=None):
    """State (d, a, b) with its sound m=3 contribution interval."""
    catalog_windows = {} if catalog_windows is None else catalog_windows
    d = _exact_int(d, "d")
    a = _exact_int(a, "a")
    b = _exact_int(b, "b")
    if d not in STATE_DEGREES:
        raise ValueError(
            f"degree {d} outside the n={N_TARGET} window "
            f"{STATE_DEGREES[0]}..{STATE_DEGREES[-1]}"
        )
    m = N_TARGET - 1 - d
    for label, order, shortfall in (("a", d, a), ("b", m, b)):
        low, high = _edge_window(edge_bounds, order)
        if not 0 <= shortfall <= high - low:
            raise ValueError(
                f"{label}={shortfall} outside the order-{order} deficiency "
                f"range 0..{high - low}"
            )
    ex = _edge_window(edge_bounds, d)[1] - a
    ey = _edge_window(edge_bounds, m)[1] - b
    return _state_from_q(
        d, a, b, ex, ey,
        q_interval(d, ex, windows, catalog_windows),
        q_interval(m, ey, windows, catalog_windows),
        budget2(d, edge_bounds),
    )


def build_states(windows, edge_bounds, catalog_windows):
    """Every state of the n=45 cone, sorted, with per-class q caching.

    A pair is omitted only when one side's edge class has no integer degree
    histogram at all. A missing R(4,5) catalog is never infeasibility: such a
    class keeps its outer interval.
    """
    intervals = {}  # explicit (order, edges) -> interval or None cache

    def resolve(order, edges):
        key = (order, edges)
        if key not in intervals:
            try:
                intervals[key] = q_interval(
                    order, edges, windows, catalog_windows
                )
            except InfeasibleEdgeClass:
                intervals[key] = None
        return intervals[key]

    states = []
    for d in STATE_DEGREES:
        m = N_TARGET - 1 - d
        x_low, x_high = _edge_window(edge_bounds, d)
        y_low, y_high = _edge_window(edge_bounds, m)
        budget = budget2(d, edge_bounds)
        x_side = [
            (a, resolve(d, x_high - a)) for a in range(x_high - x_low + 1)
        ]
        y_side = [
            (b, resolve(m, y_high - b)) for b in range(y_high - y_low + 1)
        ]
        for a, x_interval in x_side:
            if x_interval is None:
                continue
            for b, y_interval in y_side:
                if y_interval is None:
                    continue
                states.append(_state_from_q(
                    d, a, b, x_high - a, y_high - b,
                    x_interval, y_interval, budget,
                ))
    return sorted(states)


# ------------------------------------------------------ catalog provenance ---

SCHEMA_VERSION = 1
CAMPAIGN_ID = "higher_order_identity_positive_deficiency_m3"
DISPOSITION = "M3_ANALYZED"       # Task 3 owns every later disposition
OUTER_SOURCE = "outer_degree_histogram"

DATA_DIRNAME = "data"             # r55/data owns every catalog and table
TABLE_FILE = "structural_tables.json"
VALIDATION_FILE = "VALIDATION.json"
VALIDATION_EXTREME_FILE = "VALIDATION_extreme.json"
R35_FILES = tuple(f"r35_{order}.g6" for order in range(1, R35_MAX_ORDER + 1))
R45_CENSUS_FILE = "r45_24.g6"
R45_CENSUS_GRAPHS = 352366        # published complete order-24 R(4,5) census
R45_EXTREME_DIRNAME = "r45extreme"
R45_EXTREME_ORDERS = tuple(range(17, R45_MAX_ORDER))   # 17..23
R45_EXTREME_RECORDS = 36          # validated fixed-edge classes in that band
R55_REPLAY_FILE = "r55_42some.g6"
R55_REPLAY_ORDER = 42
R55_REPLAY_GRAPHS = 328           # published R(5,5,42) representatives
EXECUTED_SOURCES = (
    "check_ramsey.py", "m3_deficiency_cone.py", "subgraph_identities.py",
)

_SRC_DIR = Path(__file__).resolve().parent

# Published completeness is trusted, never re-proved here. Local hashing only
# pins *which* bytes were consumed; it says nothing about the census claims.
TRUST_ROOTS = (
    (
        "The 13 catalogs r35_1.g6..r35_13.g6 are the complete censuses of "
        "R(3,5) graphs; their completeness is trusted from the McKay data "
        "page and is not re-proved locally.",
        "https://users.cecs.anu.edu.au/~bdm/data/ramsey.html",
    ),
    (
        "Each selected order-17..23 file of data/r45extreme is the complete "
        "census of one fixed-edge R(4,5) class; completeness is trusted from "
        "the Angeltveit-McKay census of R(4,5,n,e>=e0) and its data-page "
        "archive r45extreme.tar.gz, and is not re-proved locally.",
        "https://arxiv.org/abs/2409.15709",
    ),
    (
        "r45_24.g6 is the complete census of the 352366 R(4,5,24) graphs; "
        "completeness and the count are trusted from the 2018 "
        "Angeltveit-McKay paper completing the R(4,5) extremal catalogue and "
        "its data page, and are not re-proved locally.",
        "https://arxiv.org/abs/1703.08768",
    ),
    (
        "r55_42some.g6 holds 328 published R(5,5,42) representatives; the "
        "replay proves the m=3 identity on them and their complements, not "
        "that the 328 exhaust R(5,5,42).",
        "https://users.cecs.anu.edu.au/~bdm/data/ramsey.html",
    ),
)


def _is_sha256(value):
    """True only for a lowercase 64-hex digest string."""
    return (
        isinstance(value, str)
        and len(value) == 64
        and all(char in "0123456789abcdef" for char in value)
    )


def _sha256(path):
    """Streamed SHA-256 of one consumed input file."""
    digest = hashlib.sha256()
    try:
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1 << 20), b""):
                digest.update(block)
    except OSError as exc:
        raise ValueError(f"{path.name}: unreadable input ({exc})") from None
    return digest.hexdigest()


def _relative_path(path, math_root):
    """Canonical `math/`-relative POSIX path; no absolute path is exported."""
    try:
        relative = path.resolve().relative_to(math_root)
    except (OSError, ValueError):
        raise ValueError(
            f"{path.name}: outside the math root, refusing to export a path"
        ) from None
    return relative.as_posix()


def _input_record(path, math_root, graph_count, expected_sha):
    """One provenance record; a recorded digest must match the local bytes."""
    digest = _sha256(path)
    if expected_sha is not None and digest != expected_sha:
        raise ValueError(
            f"{path.name}: local sha256 {digest} does not match the "
            f"validated {expected_sha}"
        )
    try:
        size = path.stat().st_size
    except OSError as exc:
        raise ValueError(f"{path.name}: unstatable input ({exc})") from None
    if size <= 0:
        raise ValueError(f"{path.name}: empty input")
    if graph_count is not None and _exact_int(graph_count, "graph_count") <= 0:
        raise ValueError(f"{path.name}: non-positive graph count")
    return {
        "relative_path": _relative_path(path, math_root),
        "sha256": digest,
        "bytes": size,
        "graph_count": graph_count,
    }


def _verify_unchanged_input(path, record):
    """Require the bytes hashed before analysis to remain unchanged."""
    digest = _sha256(path)
    if digest != record["sha256"]:
        raise ValueError(f"{path.name}: input changed during analysis")
    try:
        size = path.stat().st_size
    except OSError as exc:
        raise ValueError(f"{path.name}: unstatable input ({exc})") from None
    if size != record["bytes"]:
        raise ValueError(f"{path.name}: input size changed during analysis")


def _sorted_inputs(records):
    """Provenance sorted by canonical path, duplicate paths rejected."""
    ordered = sorted(records, key=lambda record: record["relative_path"])
    for earlier, later in zip(ordered, ordered[1:]):
        if earlier["relative_path"] == later["relative_path"]:
            raise ValueError(
                f"duplicate input record for {later['relative_path']}"
            )
    return ordered


# ----------------------------------------------------- validation records ---


def _record_filename(record, source):
    """Filename and parent directory name of a validation record.

    The recorded path is an absolute machine path, so only its tail is used:
    every file is reopened under the repository root derived from `__file__`.
    """
    raw = record.get("file")
    if not isinstance(raw, str) or not raw:
        raise ValueError(f"{source}: result record without a file field")
    pure = PurePosixPath(raw)
    if pure.name in ("", ".", ".."):
        raise ValueError(f"{source}: record file {raw!r} has no filename")
    return pure.name, pure.parent.name


def _require_validated(record, keys, label, source):
    """Every named predicate exactly true, and no recorded violation."""
    for key in keys:
        if record.get(key) is not True:
            raise ValueError(f"{source}: {label} does not have {key} == true")
    violations = record.get("violations")
    if not isinstance(violations, list) or violations:
        raise ValueError(f"{source}: {label} has violations {violations!r}")


def _record_order_and_count(record, label, source):
    """The single order of a record and its agreed graph count."""
    orders = record.get("n_values")
    if not isinstance(orders, dict) or len(orders) != 1:
        raise ValueError(f"{source}: {label} does not cover exactly one order")
    key, count = next(iter(orders.items()))
    if not (isinstance(key, str) and key.isascii() and key.isdigit()):
        raise ValueError(f"{source}: {label} has a non-numeric order {key!r}")
    order = int(key)
    count = _exact_int(count, f"{label} n_values[{order}]")
    graphs = _exact_int(record.get("graphs"), f"{label} graphs")
    if graphs <= 0 or graphs != count:
        raise ValueError(
            f"{source}: {label} claims {graphs} graphs but {count} at "
            f"order {order}"
        )
    return order, graphs


def _select_main_records(path):
    """The 13 R(3,5) catalogs, the order-24 census, and the replay file.

    Only fully validated single-order records carrying a recorded SHA-256 are
    accepted, so a missing or unvalidated record fails the run instead of
    silently weakening the analysis.
    """
    source = path.name
    raw = _load_json(path)
    results = raw.get("results")
    if not isinstance(results, list):
        raise ValueError(f"{source}: no results list")
    wanted = {
        name: (3, 5, order, None)
        for order, name in enumerate(R35_FILES, start=1)
    }
    wanted[R45_CENSUS_FILE] = (4, 5, R45_MAX_ORDER, R45_CENSUS_GRAPHS)
    wanted[R55_REPLAY_FILE] = (5, 5, R55_REPLAY_ORDER, R55_REPLAY_GRAPHS)
    selected = {}
    for record in results:
        if not isinstance(record, dict):
            raise ValueError(f"{source}: malformed result record")
        name, parent = _record_filename(record, source)
        if name not in wanted:
            continue
        if parent != DATA_DIRNAME:
            raise ValueError(
                f"{source}: {name} is recorded outside {DATA_DIRNAME}/"
            )
        if name in selected:
            raise ValueError(f"{source}: duplicate record for {name}")
        want_s, want_t, want_order, want_graphs = wanted[name]
        if (record.get("s"), record.get("t")) != (want_s, want_t):
            raise ValueError(
                f"{source}: {name} is not an R({want_s},{want_t}) record"
            )
        _require_validated(
            record, ("ok", "count_matches", "all_ramsey"), name, source
        )
        order, graphs = _record_order_and_count(record, name, source)
        if order != want_order:
            raise ValueError(
                f"{source}: {name} covers order {order}, expected {want_order}"
            )
        if want_graphs is not None and graphs != want_graphs:
            raise ValueError(
                f"{source}: {name} has {graphs} graphs, expected {want_graphs}"
            )
        digest = record.get("sha256")
        if not _is_sha256(digest):
            raise ValueError(f"{source}: {name} has no recorded sha256")
        selected[name] = {"order": order, "graphs": graphs, "sha256": digest}
    missing = sorted(set(wanted) - set(selected))
    if missing:
        raise ValueError(f"{source}: no validated record for {missing}")
    return selected


def _extreme_name_class(name, source):
    """The (order, edges) class encoded by an `r45<order>.<edges>.g6` name."""
    if not name.startswith("r45") or not name.endswith(".g6"):
        raise ValueError(f"{source}: {name} is not an r45 edge-class catalog")
    fields = name[3:-3].split(".")
    if len(fields) != 2 or not all(
        field.isascii() and field.isdigit() for field in fields
    ):
        raise ValueError(f"{source}: {name} does not encode order and edges")
    return int(fields[0]), int(fields[1])


def _select_extreme_records(path):
    """Every validated fixed-edge R(4,5) class at orders 17..23.

    This report has no `count_matches` or SHA field, so neither is required;
    the local digest is computed and persisted instead. Filename, metadata
    order, and metadata edge class must agree, and the band must be complete.
    """
    source = path.name
    raw = _load_json(path)
    results = raw.get("results")
    if not isinstance(results, list):
        raise ValueError(f"{source}: no results list")
    selected = {}
    for record in results:
        if not isinstance(record, dict):
            raise ValueError(f"{source}: malformed result record")
        name, parent = _record_filename(record, source)
        order, graphs = _record_order_and_count(record, name, source)
        if order not in R45_EXTREME_ORDERS:
            continue
        if parent != R45_EXTREME_DIRNAME:
            raise ValueError(
                f"{source}: target-order record {name} is outside "
                f"{R45_EXTREME_DIRNAME}/"
            )
        if (record.get("s"), record.get("t")) != (4, 5):
            raise ValueError(f"{source}: {name} is not an R(4,5) record")
        _require_validated(record, ("ok", "all_ramsey"), name, source)
        edge_min, edge_max = _window(record, "edge_min", "edge_max", name)
        if edge_min != edge_max:
            raise ValueError(
                f"{source}: {name} spans edges {edge_min}..{edge_max}, not "
                "one fixed class"
            )
        if _extreme_name_class(name, source) != (order, edge_min):
            raise ValueError(
                f"{source}: {name} disagrees with its metadata "
                f"(order {order}, edges {edge_min})"
            )
        if name in selected:
            raise ValueError(f"{source}: duplicate record for {name}")
        selected[name] = {"order": order, "edges": edge_min, "graphs": graphs}
    if len(selected) != R45_EXTREME_RECORDS:
        raise ValueError(
            f"{source}: {len(selected)} validated order-"
            f"{R45_EXTREME_ORDERS[0]}..{R45_EXTREME_ORDERS[-1]} records, "
            f"expected {R45_EXTREME_RECORDS}"
        )
    covered = {info["order"] for info in selected.values()}
    missing = [order for order in R45_EXTREME_ORDERS if order not in covered]
    if missing:
        raise ValueError(f"{source}: no edge class for orders {missing}")
    return selected


# -------------------------------------------------------- catalog streams ---


def _normalized_catalog_record(record):
    """Validate one `(path, order, edges_or_none, count)` catalog record."""
    if type(record) is not tuple:
        raise ValueError(f"catalog record must be a 4-tuple, got {record!r}")
    if len(record) != 4:
        raise ValueError(
            f"catalog record must have 4 fields, got {len(record)}"
        )
    path, order, edges, count = record
    if not isinstance(path, (str, Path)):
        raise ValueError(f"catalog path must be a path, got {path!r}")
    path = Path(path)
    order = _exact_int(order, f"{path.name} expected order")
    if order <= 0:
        raise ValueError(f"{path.name}: expected order must be positive")
    if edges is not None:
        edges = _exact_int(edges, f"{path.name} expected edges")
        if edges < 0:
            raise ValueError(f"{path.name}: expected edges must be >= 0")
    count = _exact_int(count, f"{path.name} expected graph count")
    if count <= 0:
        raise ValueError(
            f"{path.name}: expected graph count must be positive, got {count}"
        )
    return path, order, edges, count


def _stream_q_counts(path, expected_order, expected_edges, expected_count):
    """Exact q frequencies per (order, edges) for one graph6 catalog.

    One pass, no parsed graph retained: memory stays proportional to the tiny
    frequency maps. `motif3` is never called; the only cubic-ish work is the
    single `triangle_count` each graph needs.
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
                wedges = 0
                degree_sum = 0
                for row in adj:
                    degree = popcount(row)
                    degree_sum += degree
                    wedges += degree * (degree - 1) // 2
                edges = degree_sum // 2
                if expected_edges is not None and edges != expected_edges:
                    raise ValueError(
                        f"{path.name} line {number}: {edges} edges, expected "
                        f"{expected_edges}"
                    )
                q = wedges - triangle_count(adj)
                frequency = counts.get((order, edges))
                if frequency is None:
                    frequency = counts[(order, edges)] = {}
                frequency[q] = frequency.get(q, 0) + 1
        except (OSError, UnicodeError) as exc:
            raise ValueError(
                f"{path.name}: unreadable catalog ({exc})"
            ) from None
    if not graphs:
        raise ValueError(f"{path.name}: empty catalog")
    if graphs != expected_count:
        raise ValueError(
            f"{path.name}: {graphs} graphs, expected {expected_count}"
        )
    return graphs, {
        key: {q: frequency[q] for q in sorted(frequency)}
        for key, frequency in sorted(counts.items())
    }


def _catalog_sweep(catalog_records):
    """Stream every catalog once: class histograms, line counts, owners."""
    records = list(catalog_records)
    if not records:
        raise ValueError("no catalog records to stream")
    histograms = {}
    graph_counts = {}
    owners = {}
    for record in records:
        path, order, edges, count = _normalized_catalog_record(record)
        location = path.as_posix()
        if location in graph_counts:
            raise ValueError(f"catalog {path.name} is listed twice")
        graphs, produced = _stream_q_counts(path, order, edges, count)
        graph_counts[location] = graphs
        for key, frequency in produced.items():
            owner = owners.get(key)
            if owner is not None:
                raise ValueError(
                    f"class (order={key[0]}, edges={key[1]}) is produced by "
                    f"both {owner.name} and {path.name}"
                )
            owners[key] = path
            histograms[key] = frequency
    return histograms, graph_counts, owners


def catalog_q_histograms(catalog_records):
    """Exact q frequency map per (order, edges) over streamed catalogs.

    `catalog_records` is an iterable of normalized 4-tuples
    `(path, expected_order, expected_edges_or_none, expected_graph_count)`.
    Every non-empty graph6 line is parsed, its order checked against
    `expected_order`, its exact edge count checked against `expected_edges`
    when the class is fixed, and

        q = sum_v binom(deg v, 2) - triangle_count

    accumulated into the frequency map of its own `(order, edges)` class. The
    final line count must equal `expected_graph_count`.

    Fail-closed: a malformed record, a non-positive or non-integer count, an
    empty file, a wrong order, a wrong fixed edge count, a wrong total, or an
    `(order, edges)` class produced by two different source files raises
    `ValueError`. Keys are inserted in sorted order, per class and per q.
    """
    return _catalog_sweep(catalog_records)[0]


def _count_graph6_lines(path):
    """Non-empty graph6 line count of one small catalog."""
    try:
        raw = path.read_bytes()
    except OSError as exc:
        raise ValueError(f"{path.name}: unreadable catalog ({exc})") from None
    count = sum(1 for line in raw.splitlines() if line.strip())
    if not count:
        raise ValueError(f"{path.name}: empty catalog")
    return count


def _replay_published(path, expected_order, expected_graphs):
    """Zero scaled m=3 residual on every published graph and its complement."""
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
                    residual = m3_residual_scaled(graph)
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
            f"{path.name}: {graphs} graphs, expected {expected_graphs}"
        )
    if zeros != 2 * expected_graphs:
        raise ValueError(
            f"{path.name}: {zeros} zero residuals, expected "
            f"{2 * expected_graphs}"
        )
    return graphs, zeros


# --------------------------------------------------------------- analysis ---


def _state_record(state, edge_bounds, class_source):
    """One canonical state record, with per-side q interval provenance."""
    m = N_TARGET - 1 - state.d
    x_edges = _edge_window(edge_bounds, state.d)[1] - state.a
    y_edges = _edge_window(edge_bounds, m)[1] - state.b
    x_source = class_source.get((state.d, x_edges), OUTER_SOURCE)
    y_source = class_source.get((m, y_edges), OUTER_SOURCE)
    return {
        "d": state.d,
        "a": state.a,
        "b": state.b,
        "deficiency": state.deficiency,
        "excess_balance": state.excess_balance,
        "g_lo": state.g_lo,
        "g_hi": state.g_hi,
        "x_interval_source": x_source,
        "y_interval_source": y_source,
    }


def _catalog_windows(histograms):
    """Exact q window of every nonempty class, in sorted class order."""
    windows = {}
    for key, frequency in sorted(histograms.items()):
        if not frequency:
            raise ValueError(
                f"class (order={key[0]}, edges={key[1]}) has no q values"
            )
        windows[key] = (min(frequency), max(frequency))
    return windows


def _document(
    inputs, windows, histograms, class_source, states, edge_bounds,
    replay_graphs, replay_zeros,
):
    """The canonical document; every list is deterministically ordered."""
    catalog_windows = _catalog_windows(histograms)
    return {
        "schema_version": SCHEMA_VERSION,
        "campaign_id": CAMPAIGN_ID,
        "inputs": _sorted_inputs(inputs),
        "trust_roots": [
            {"statement": statement, "url": url}
            for statement, url in sorted(TRUST_ROOTS)
        ],
        "replay": {
            "published_graphs": replay_graphs,
            "complements": replay_graphs,
            "residual_zero": replay_zeros,
        },
        "r35_edge_windows": [
            {
                "order": order,
                "edge_min": windows[order][0],
                "edge_max": windows[order][1],
            }
            for order in sorted(windows)
        ],
        "catalog_q_histograms": [
            {
                "order": key[0],
                "edges": key[1],
                "source": class_source[key],
                "graph_count": sum(histograms[key].values()),
                "q_counts": [
                    {"q": q, "count": histograms[key][q]}
                    for q in sorted(histograms[key])
                ],
            }
            for key in sorted(histograms)
        ],
        "catalog_windows": [
            {
                "order": key[0],
                "edges": key[1],
                "q_min": catalog_windows[key][0],
                "q_max": catalog_windows[key][1],
                "source": class_source[key],
            }
            for key in sorted(catalog_windows)
        ],
        "states": [
            _state_record(state, edge_bounds, class_source)
            for state in states
        ],
        "searches": [],
        "disposition": DISPOSITION,
    }


def _analyze(src_dir, echo):
    """The whole fail-closed analysis; returns the canonical document."""
    data_dir = src_dir.parent / DATA_DIRNAME
    math_root = src_dir.parents[1]
    validation = data_dir / VALIDATION_FILE
    extreme_validation = data_dir / VALIDATION_EXTREME_FILE
    tables = data_dir / TABLE_FILE

    edge_bounds = load_edge_bounds(tables)
    windows = r35_edge_windows(validation)
    main_records = _select_main_records(validation)
    extreme_records = _select_extreme_records(extreme_validation)

    # Bind every consumed path before any expensive computation. Re-hashing
    # the same records after the sweeps closes the provenance TOCTOU window:
    # the JSON can describe only the exact bytes used throughout this run.
    inputs = []
    preflight = {}

    def add_input(path, graph_count, expected_sha=None):
        record = _input_record(
            path, math_root, graph_count, expected_sha
        )
        inputs.append(record)
        preflight[path] = record

    for name in R35_FILES:
        record = main_records[name]
        path = data_dir / name
        graphs = _count_graph6_lines(path)
        if graphs != record["graphs"]:
            raise ValueError(
                f"{name}: {graphs} graphs, expected {record['graphs']}"
            )
        add_input(path, graphs, record["sha256"])
    for name in sorted(extreme_records):
        record = extreme_records[name]
        add_input(
            data_dir / R45_EXTREME_DIRNAME / name, record["graphs"]
        )
    census = main_records[R45_CENSUS_FILE]
    census_path = data_dir / R45_CENSUS_FILE
    add_input(census_path, census["graphs"], census["sha256"])
    replay_path = data_dir / R55_REPLAY_FILE
    add_input(
        replay_path, main_records[R55_REPLAY_FILE]["graphs"],
        main_records[R55_REPLAY_FILE]["sha256"],
    )
    for path in (validation, extreme_validation, tables):
        add_input(path, None)
    for name in EXECUTED_SOURCES:
        add_input(src_dir / name, None)

    inputs = _sorted_inputs(inputs)

    replay_path = data_dir / R55_REPLAY_FILE
    replay_graphs, replay_zeros = _replay_published(
        replay_path, R55_REPLAY_ORDER, R55_REPLAY_GRAPHS
    )
    echo(f"M3 REPLAY: {replay_zeros}/{2 * replay_graphs} residual zero")

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
    census = main_records[R45_CENSUS_FILE]
    census_path = data_dir / R45_CENSUS_FILE
    # The complete census spans every extremal class, so its edges are free.
    catalog_records.append(
        (census_path, census["order"], None, census["graphs"])
    )
    histograms, graph_counts, owners = _catalog_sweep(catalog_records)

    outer = {}          # one outer interval per (order, edges), cached
    violations = 0
    for key in sorted(histograms):
        if key not in outer:
            try:
                outer[key] = q_outer_interval(key[0], key[1], windows)
            except InfeasibleEdgeClass:
                raise ValueError(
                    f"catalog class (order={key[0]}, edges={key[1]}) has no "
                    "integer degree histogram"
                ) from None
        low, high = outer[key]
        violations += sum(1 for q in histograms[key] if not low <= q <= high)
    if violations:
        raise ValueError(
            f"{violations} exact catalog q values fall outside the outer "
            "degree-histogram envelope"
        )
    echo(f"M3 ENVELOPES: catalog violations={violations}")

    catalog_windows = _catalog_windows(histograms)
    states = build_states(windows, edge_bounds, catalog_windows)
    if not states:
        raise ValueError("the n=45 state cone is empty")
    echo(f"M3 STATES: {len(states)} states")

    class_source = {
        key: _relative_path(path, math_root) for key, path in owners.items()
    }
    for path, record in preflight.items():
        _verify_unchanged_input(path, record)

    return _document(
        inputs, windows, histograms, class_source, states, edge_bounds,
        replay_graphs, replay_zeros,
    )


# ------------------------------------------------- serialization and CLI ---

_TOP_LEVEL_KEYS = (
    "campaign_id", "catalog_q_histograms", "catalog_windows", "disposition",
    "inputs", "r35_edge_windows", "replay", "schema_version", "searches",
    "states", "trust_roots",
)
_INPUT_KEYS = ("bytes", "graph_count", "relative_path", "sha256")
_TRUST_KEYS = ("statement", "url")
_REPLAY_KEYS = ("complements", "published_graphs", "residual_zero")
_R35_WINDOW_KEYS = ("edge_max", "edge_min", "order")
_HISTOGRAM_KEYS = ("edges", "graph_count", "order", "q_counts", "source")
_CATALOG_WINDOW_KEYS = ("edges", "order", "q_max", "q_min", "source")
_Q_COUNT_KEYS = ("count", "q")
_STATE_KEYS = (
    "a", "b", "d", "deficiency", "excess_balance", "g_hi", "g_lo",
    "x_interval_source", "y_interval_source",
)
_EXPECTED_INPUTS = (
    len(R35_FILES)              # the complete R(3,5) catalogs
    + R45_EXTREME_RECORDS       # the selected fixed-edge classes
    + 2                         # the order-24 census and the replay file
    + 3                         # both validation reports and the edge tables
    + len(EXECUTED_SOURCES)     # the executed Python sources
)


def _require_keys(record, keys, label):
    """Exactly the schema keys, nothing missing and nothing extra."""
    if not isinstance(record, dict):
        raise ValueError(f"{label} is not an object")
    if tuple(sorted(record)) != tuple(keys):
        raise ValueError(
            f"{label} has keys {sorted(record)}, expected {list(keys)}"
        )
    return record


def _iter_strings(value):
    """Every string in the document, keys included."""
    if isinstance(value, str):
        yield value
    elif isinstance(value, dict):
        for key, item in value.items():
            yield key
            for text in _iter_strings(item):
                yield text
    elif isinstance(value, list):
        for item in value:
            for text in _iter_strings(item):
                yield text


def _validate_document(parsed, document):
    """Schema and count invariants, re-checked on the parsed-back file."""
    if parsed != document:
        raise ValueError("the serialized analysis does not parse back exactly")
    _require_keys(parsed, _TOP_LEVEL_KEYS, "analysis document")
    if parsed["schema_version"] != SCHEMA_VERSION:
        raise ValueError(f"schema_version {parsed['schema_version']!r}")
    if parsed["campaign_id"] != CAMPAIGN_ID:
        raise ValueError(f"campaign_id {parsed['campaign_id']!r}")
    if parsed["disposition"] != DISPOSITION:
        raise ValueError(f"disposition {parsed['disposition']!r}")
    if parsed["searches"] != []:
        raise ValueError("searches must stay empty until Task 3")
    for text in _iter_strings(parsed):
        if text.startswith("/"):
            raise ValueError(f"absolute path {text!r} in the analysis")

    if _require_keys(parsed["replay"], _REPLAY_KEYS, "replay") != {
        "published_graphs": R55_REPLAY_GRAPHS,
        "complements": R55_REPLAY_GRAPHS,
        "residual_zero": 2 * R55_REPLAY_GRAPHS,
    }:
        raise ValueError("replay is not the full published-plus-dual sweep")

    inputs = parsed["inputs"]
    if not isinstance(inputs, list):
        raise ValueError("inputs is not a list")
    if len(inputs) != _EXPECTED_INPUTS:
        raise ValueError(
            f"{len(inputs)} input records, expected {_EXPECTED_INPUTS}"
        )
    paths = []
    for record in inputs:
        _require_keys(record, _INPUT_KEYS, "input record")
        label = record["relative_path"]
        if not isinstance(label, str) or not label:
            raise ValueError("input record without a relative path")
        if not _is_sha256(record["sha256"]):
            raise ValueError(f"{label}: malformed sha256")
        if _exact_int(record["bytes"], f"{label} bytes") <= 0:
            raise ValueError(f"{label}: empty input")
        count = record["graph_count"]
        if count is not None and _exact_int(count, f"{label} count") <= 0:
            raise ValueError(f"{label}: non-positive graph count")
        paths.append(label)
    if paths != sorted(set(paths)):
        raise ValueError("input records are not sorted and unique")

    roots = parsed["trust_roots"]
    if not isinstance(roots, list) or not roots:
        raise ValueError("trust_roots must be a non-empty list")
    for record in roots:
        _require_keys(record, _TRUST_KEYS, "trust root")
        if not record["statement"]:
            raise ValueError("trust root without a statement")
        if not str(record["url"]).startswith("https://"):
            raise ValueError(f"trust root url {record['url']!r}")
    if roots != sorted(roots, key=lambda r: (r["statement"], r["url"])):
        raise ValueError("trust_roots are not sorted")

    r35_windows = parsed["r35_edge_windows"]
    if not isinstance(r35_windows, list):
        raise ValueError("r35_edge_windows is not a list")
    orders = []
    for record in r35_windows:
        _require_keys(record, _R35_WINDOW_KEYS, "R(3,5) window")
        if record["edge_min"] > record["edge_max"]:
            raise ValueError(f"empty R(3,5) window at order {record['order']}")
        orders.append(record["order"])
    if orders != list(range(R35_MAX_ORDER + 1)):
        raise ValueError(f"R(3,5) windows cover {orders}")

    histograms = parsed["catalog_q_histograms"]
    if not isinstance(histograms, list):
        raise ValueError("catalog_q_histograms is not a list")
    if len(histograms) <= R45_EXTREME_RECORDS:
        raise ValueError(
            f"{len(histograms)} catalog classes, expected more than "
            f"{R45_EXTREME_RECORDS}"
        )
    classes = []
    sources = set()
    observed = {}
    for record in histograms:
        _require_keys(record, _HISTOGRAM_KEYS, "catalog histogram")
        key = (record["order"], record["edges"])
        classes.append(key)
        source = record["source"]
        if not isinstance(source, str) or not source.endswith(".g6"):
            raise ValueError(f"class {key} has no graph6 source")
        sources.add(source)
        entries = record["q_counts"]
        if not isinstance(entries, list) or not entries:
            raise ValueError(f"class {key} has an empty q histogram")
        values = []
        total = 0
        for entry in entries:
            _require_keys(entry, _Q_COUNT_KEYS, "q count")
            frequency = _exact_int(entry["count"], f"class {key} count")
            if frequency <= 0:
                raise ValueError(f"class {key} has a non-positive frequency")
            values.append(_exact_int(entry["q"], f"class {key} q"))
            total += frequency
        if values != sorted(set(values)):
            raise ValueError(f"class {key} has unsorted or repeated q values")
        graph_count = _exact_int(record["graph_count"], f"class {key} graphs")
        if total != graph_count or graph_count <= 0:
            raise ValueError(
                f"class {key} sums to {total}, not its {graph_count} graphs"
            )
        observed[key] = (values[0], values[-1], source)
    if classes != sorted(set(classes)):
        raise ValueError("catalog histograms are not sorted and unique")

    exact_windows = parsed["catalog_windows"]
    if not isinstance(exact_windows, list):
        raise ValueError("catalog_windows is not a list")
    if len(exact_windows) != len(classes):
        raise ValueError(
            f"{len(exact_windows)} catalog windows for {len(classes)} classes"
        )
    seen = []
    for record in exact_windows:
        _require_keys(record, _CATALOG_WINDOW_KEYS, "catalog window")
        key = (record["order"], record["edges"])
        if key not in observed:
            raise ValueError(f"catalog window {key} has no histogram")
        if (record["q_min"], record["q_max"], record["source"]) != observed[
            key
        ]:
            raise ValueError(
                f"catalog window {key} disagrees with its own histogram"
            )
        seen.append(key)
    if seen != sorted(set(seen)):
        raise ValueError("catalog windows are not sorted and unique")

    states = parsed["states"]
    if not isinstance(states, list) or not states:
        raise ValueError("the analysis records no states")
    cells = []
    for record in states:
        _require_keys(record, _STATE_KEYS, "state")
        if record["deficiency"] != record["a"] + record["b"]:
            raise ValueError("state deficiency disagrees with a + b")
        if record["g_lo"] > record["g_hi"]:
            raise ValueError("state has an empty g interval")
        for field in ("x_interval_source", "y_interval_source"):
            value = record[field]
            if value != OUTER_SOURCE and value not in sources:
                raise ValueError(f"state {field} {value!r} is not a catalog")
        cells.append((record["d"], record["a"], record["b"]))
    if cells != sorted(set(cells)):
        raise ValueError("states are not sorted and unique")


def _write_document(target, document):
    """Serialize deterministically, parse back, then atomically replace."""
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
        os.chmod(name, 0o644)   # mkstemp is 0600; artifacts are world-readable
        _validate_document(_load_json(temporary), document)
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
    """Calibrate the cone against the catalogs and write the analysis JSON."""
    parser = argparse.ArgumentParser(
        description=(
            "Calibrate the m=3 deficiency cone against the published R(4,5) "
            "catalogs and write the canonical analysis JSON."
        )
    )
    parser.add_argument(
        "--output", required=True,
        help="canonical analysis JSON, replaced atomically after every check",
    )
    args = parser.parse_args(argv)
    try:
        document = _analyze(_SRC_DIR, _echo)
        _write_document(Path(args.output), document)
    except (ValueError, OSError) as exc:
        print(f"m3 analysis failed: {exc}", file=sys.stderr)
        return 1
    print(f"WROTE {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
