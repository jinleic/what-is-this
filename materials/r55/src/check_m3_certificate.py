#!/usr/bin/env python3
"""Independent checker for the m=3 higher-identity certificate artifact.

Role
----
`search_m3_cuts.py` writes `data/higher_identity_m3.json` as the single
certificate artifact of the m=3 campaign.  This module re-verifies that
artifact from scratch, sharing no code with the producer chain: it imports
neither `subgraph_identities`, `m3_deficiency_cone`, nor `search_m3_cuts`.
Every kernel primitive, selection rule, envelope, state, count, hash,
certificate, witness, and disposition below is recomputed here from the
frozen constants and the recorded input bytes, and any disagreement is a
`CheckViolation` naming the offending field.

What is checked, in fail-fast stage order
-----------------------------------------
1.  Schema preflight: exact key sets everywhere, canonical `r55/...`
    relative paths, no absolute paths anywhere in the document, sorted
    unique lists, per-section self-consistency (histogram sums, catalog
    windows against their own histograms, state invariants).
2.  Hash preflight: every recorded input exists under the root and its
    streamed SHA-256 equals the recorded value.
3.  Replay: the 328 published R(5,5,42) graphs and their complements are
    streamed and the checker's OWN scaled m=3 residual kernel must return
    zero on every one; the replay section must record exactly what the
    file yields.
4.  Search evidence: the four search records are verified against the
    artifact's OWN stored states, so a corrupted coefficient, bound,
    weight, or status is caught before the expensive catalog rescan.
5.  Catalog rescan: validation records are re-selected under the frozen
    rules (r44 decoys skipped, exactly 36 extreme classes at orders
    17..23, r35 window coverage 1..13), every recorded byte count equals
    the file size on disk, every catalog graph is re-parsed and
    re-counted, the replay file's line count is tied back to its own
    validation record, the observed per-class q histograms, owners,
    windows, and input set must equal the recorded ones, the census must
    own the whole observed order-24 band, and every observed q must lie
    inside the checker's own outer degree-histogram envelope.
6.  State rebuild: windows, edge bounds, and the cone are rebuilt by the
    sole g formula and must equal every recorded state field, including
    the per-side interval provenance strings.
7.  Terminal mapping: the disposition must equal the verdict implied by
    the three route statuses.
8.  Final re-hash: every recorded input is hashed again after the whole
    analysis, so a mid-run mutation of any file is still caught (the
    producer closes the same window after its own long run).

CLI
---
`check_m3_certificate.py [--root DIR] ANALYSIS` exits 0 and prints
`M3 EVIDENCE VERIFIED: <disposition>` on the last stdout line, or exits 1
with `M3 EVIDENCE FAILED: <field>: <counterexample>` on stderr.  `--root`
substitutes the math repository root that every recorded `r55/...` path
resolves against; the default is the checker's own `parents[2]`.

Performance: the rescan streams 8,500,211 catalog graphs once, single
process, parsing each line, popcounting 24 rows, and counting triangles
with pure bit operations; no graph is retained.
"""

import argparse
import hashlib
import json
import os
import sys
from fractions import Fraction
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_ramsey import parse_graph6_line, popcount  # noqa: E402


# ------------------------------------------------------- frozen constants ---

N_TARGET = 45        # order of the hypothetical R(5,5,45) graph
R44 = 18             # R(4,4): lower bound of an R(4,5) degree window
R35_MAX_ORDER = 13   # R(3,5) - 1: upper bound of an R(4,5) degree window
R45_MAX_ORDER = 24   # R(4,5) - 1: maximum order of an R(4,5) graph
STATE_DEGREES = tuple(
    range(N_TARGET - 1 - R45_MAX_ORDER, R45_MAX_ORDER + 1)
)
VERTICES = Fraction(N_TARGET)   # every vertex occupies exactly one state

SCHEMA_VERSION = 1
CAMPAIGN_ID = "higher_order_identity_positive_deficiency_m3"

DATA_DIRNAME = "data"           # r55/data owns every catalog and table
TABLE_FILE = "structural_tables.json"
VALIDATION_FILE = "VALIDATION.json"
VALIDATION_EXTREME_FILE = "VALIDATION_extreme.json"
R35_FILES = tuple(
    f"r35_{order}.g6" for order in range(1, R35_MAX_ORDER + 1)
)
R45_CENSUS_FILE = "r45_24.g6"
CENSUS_ORDER = R45_MAX_ORDER
CENSUS_EDGE_BAND = (116, 132)   # published extremal edge band at order 24
R45_EXTREME_DIRNAME = "r45extreme"
R45_EXTREME_ORDERS = tuple(range(17, R45_MAX_ORDER))   # 17..23
R45_EXTREME_RECORDS = 36        # validated fixed-edge classes in that band
R55_REPLAY_FILE = "r55_42some.g6"
R55_REPLAY_ORDER = 42
# The executed sources pinned by the two producing tasks; the artifact's
# input list must record exactly these four Python files.
RECORDED_SOURCES = (
    "check_ramsey.py", "m3_deficiency_cone.py", "subgraph_identities.py",
    "search_m3_cuts.py",
)
OUTER_SOURCE = "outer_degree_histogram"

TRIANGLE_COEFFICIENT = 6        # the scaled m=3 identity's triangle weight

ACCEPTED = "ACCEPTED_EXACT_CUT"
REJECTED = "REJECTED_BY_EXACT_WITNESS"
UNRESOLVED = "CERTIFICATION_UNRESOLVED"
UNAVAILABLE = "UNAVAILABLE_NO_COVER_CERTIFICATE"
ROUTE_STATUSES = (ACCEPTED, REJECTED, UNRESOLVED)
ROUTE_IDS = ("total_deficiency", "degree20_count",
             "deficiency_ge8_count")
UNAVAILABLE_ROUTE = "required_local_family"

DISPOSITION_ACCEPTED = "M3_ACCEPTED_CUT"
DISPOSITION_NO_CUT = "M3_NO_CUT_IN_FROZEN_CONE"
DISPOSITION_UNRESOLVED = "M3_CERTIFICATION_UNRESOLVED"

UNAVAILABLE_RECORD = {
    "objective_id": UNAVAILABLE_ROUTE,
    "certificate": None,
    "exact_upper_bound": None,
    "primal_witness": None,
    "exact_witness_value": None,
    "route_status": UNAVAILABLE,
}

_DEFAULT_ROOT = Path(__file__).resolve().parents[2]
_STREAM_BLOCK = 1 << 20


class CheckViolation(ValueError):
    """One artifact field disagrees with the independently derived truth."""


class InfeasibleEdgeClass(ValueError):
    """No integer degree histogram survives for an (order, edges) class."""


# ------------------------------------------------------ kernel primitives ---
# Re-derived here from scratch; the producer chain is never imported.

def complement(adj):
    """Complement graph, same vertex order, no loops."""
    n = len(adj)
    full = (1 << n) - 1
    return [full & ~adj[v] & ~(1 << v) for v in range(n)]


def induced(adj, vertices):
    """Subgraph induced on `vertices`, relabeled in the given order."""
    out = [0] * len(vertices)
    for i, u in enumerate(vertices):
        row = adj[u]
        for j in range(i + 1, len(vertices)):
            if (row >> vertices[j]) & 1:
                out[i] |= 1 << j
                out[j] |= 1 << i
    return out


def edge_count(adj):
    """Number of edges, half the handshake sum."""
    return sum(popcount(row) for row in adj) // 2


def triangle_count(adj):
    """Triangles by summing common neighbours over every edge, then //3."""
    total = 0
    for u in range(len(adj)):
        row = adj[u]
        later = row & ~((1 << (u + 1)) - 1)   # neighbours v with v > u
        while later:
            low = later & -later
            v = low.bit_length() - 1
            total += popcount(row & adj[v])
            later ^= low
    return total // 3


def induced_p3_count(adj):
    """Induced three-vertex paths: wedges minus three per triangle."""
    wedges = 0
    for row in adj:
        degree = popcount(row)
        wedges += degree * (degree - 1) // 2
    return wedges - 3 * triangle_count(adj)


def q_of_adj(adj):
    """q = sum_x binom(deg x, 2) - t, one popcount pass plus triangles."""
    wedges = 0
    for row in adj:
        degree = popcount(row)
        wedges += degree * (degree - 1) // 2
    return wedges - triangle_count(adj)


def kernel_residual_scaled(adj, tri_coeff=TRIANGLE_COEFFICIENT):
    """Total scaled m=3 residual of the checker's own identity kernel.

    Per vertex v of degree d with neighborhood X and dual G_v^-, the
    identity reads 3 t(G_v^-) = (n + 3 - 3 d) e(X) + 3 p3(X) + 6 t(X);
    the residual is the difference of the two sides, summed over all
    vertices.  It is exactly zero for EVERY graph; `tri_coeff` exists so
    a mutated coefficient can be shown to break the identity.
    """
    n = len(adj)
    full = (1 << n) - 1
    total = 0
    for v in range(n):
        row = adj[v]
        neighborhood = induced(
            adj, [u for u in range(n) if (row >> u) & 1]
        )
        dual = induced(
            adj, [
                u for u in range(n)
                if not ((row >> u) & 1) and u != v
            ]
        )
        d = len(neighborhood)
        if d + len(dual) != n - 1:
            raise CheckViolation(
                f"kernel: vertex {v} of an order-{n} graph splits "
                f"{d} + {len(dual)} != {n - 1}"
            )
        total += 3 * triangle_count(dual) - (
            (n + 3 - 3 * d) * edge_count(neighborhood)
            + 3 * induced_p3_count(neighborhood)
            + tri_coeff * triangle_count(neighborhood)
        )
    return total


# --------------------------------------------------------- value guards ---

def _exact_int(value, label):
    """Reject anything that is not already an exact Python int."""
    if type(value) is not int:
        raise CheckViolation(f"{label}: {value!r} is not an exact integer")
    return value


def _canonical_fraction(text, label):
    """The value of an already reduced canonical fraction string.

    One round trip rejects padded or signed text, decimal strings,
    unreduced fractions, and `n/1` spellings in a single stroke.
    """
    if type(text) is not str:
        raise CheckViolation(f"{label}: {text!r} must be a string")
    try:
        value = Fraction(text)
    except (ValueError, ZeroDivisionError):
        raise CheckViolation(
            f"{label}: {text!r} is not a rational"
        ) from None
    if str(value) != text:
        raise CheckViolation(
            f"{label}: {text!r} is not canonical, want {value}"
        )
    return value


def _require_keys(record, keys, label):
    """Exactly the schema keys, nothing missing and nothing extra."""
    if not isinstance(record, dict):
        raise CheckViolation(f"{label}: not an object")
    if tuple(sorted(record)) != tuple(keys):
        raise CheckViolation(
            f"{label}: has keys {sorted(record)}, expected {list(keys)}"
        )
    return record


def _is_sha256(value):
    """True only for a lowercase 64-hex digest string."""
    return (
        isinstance(value, str)
        and len(value) == 64
        and all(char in "0123456789abcdef" for char in value)
    )


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


def _load_json(path, label):
    """One parsed JSON object, or a field-naming violation."""
    try:
        raw = json.loads(path.read_text(encoding="ascii"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise CheckViolation(
            f"{label}: unreadable JSON ({exc})"
        ) from None
    if not isinstance(raw, dict):
        raise CheckViolation(f"{label}: top-level JSON object expected")
    return raw


def _sha256_of(path):
    """Streamed SHA-256 of one recorded input file."""
    digest = hashlib.sha256()
    try:
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(_STREAM_BLOCK), b""):
                digest.update(block)
    except OSError as exc:
        raise CheckViolation(
            f"sha256: {path.name}: unreadable input ({exc})"
        ) from None
    return digest.hexdigest()


# ---------------------------------------------------- validation selection ---

def _record_filename(record, source):
    """Filename and parent directory name of a validation record.

    The recorded `file` path is an absolute machine path of the producing
    run, so only its last two components are ever consumed.
    """
    raw = record.get("file")
    if not isinstance(raw, str) or not raw:
        raise CheckViolation(f"{source}: result record without a file field")
    pure = PurePosixPath(raw)
    if pure.name in ("", ".", ".."):
        raise CheckViolation(f"{source}: record file {raw!r} has no filename")
    return pure.name, pure.parent.name


def _require_validated(record, keys, label, source):
    """Every named predicate exactly true, and no recorded violation."""
    for key in keys:
        if record.get(key) is not True:
            raise CheckViolation(
                f"{source}: {label} does not have {key} == true"
            )
    violations = record.get("violations")
    if not isinstance(violations, list) or violations:
        raise CheckViolation(
            f"{source}: {label} has violations {violations!r}"
        )


def _record_order_and_count(record, label, source):
    """The single order of a record and its agreed graph count."""
    orders = record.get("n_values")
    if not isinstance(orders, dict) or len(orders) != 1:
        raise CheckViolation(
            f"{source}: {label} does not cover exactly one order"
        )
    key, count = next(iter(orders.items()))
    if not (isinstance(key, str) and key.isascii() and key.isdigit()):
        raise CheckViolation(
            f"{source}: {label} has a non-numeric order {key!r}"
        )
    order = int(key)
    count = _exact_int(count, f"{label} n_values[{order}]")
    graphs = _exact_int(record.get("graphs"), f"{label} graphs")
    if graphs <= 0 or graphs != count:
        raise CheckViolation(
            f"{source}: {label} claims {graphs} graphs but {count} at "
            f"order {order}"
        )
    return order, graphs


def select_main_records(path):
    """The validated records of the catalogs, census, and replay file.

    Wanted files: the 13 R(3,5) catalogs, the order-24 census, and the
    R(5,5,42) replay file, all under `r55/data/`.  r44_* records in the
    same report are (4,4) decoys and are skipped.  Graph counts are read
    from the report, never pinned to published values: the rescan must
    reproduce them from the bytes.
    """
    source = path.name
    raw = _load_json(path, source)
    results = raw.get("results")
    if not isinstance(results, list):
        raise CheckViolation(f"{source}: no results list")
    wanted = {
        name: (3, 5, order) for order, name in enumerate(R35_FILES, start=1)
    }
    wanted[R45_CENSUS_FILE] = (4, 5, CENSUS_ORDER)
    wanted[R55_REPLAY_FILE] = (5, 5, R55_REPLAY_ORDER)
    selected = {}
    for record in results:
        if not isinstance(record, dict):
            raise CheckViolation(f"{source}: malformed result record")
        name, parent = _record_filename(record, source)
        if name not in wanted:
            continue            # the (4,4) decoys and anything unknown
        if parent != DATA_DIRNAME:
            raise CheckViolation(
                f"{source}: {name} is recorded outside {DATA_DIRNAME}/"
            )
        if name in selected:
            raise CheckViolation(
                f"duplicate: {source}: duplicate record for {name}"
            )
        want_s, want_t, want_order = wanted[name]
        if (record.get("s"), record.get("t")) != (want_s, want_t):
            raise CheckViolation(
                f"{source}: {name} is not an R({want_s},{want_t}) record"
            )
        _require_validated(
            record, ("ok", "count_matches", "all_ramsey"), name, source
        )
        order, graphs = _record_order_and_count(record, name, source)
        if order != want_order:
            raise CheckViolation(
                f"{source}: {name} covers order {order}, expected "
                f"{want_order}"
            )
        digest = record.get("sha256")
        if not _is_sha256(digest):
            raise CheckViolation(f"{source}: {name} has no recorded sha256")
        selected[name] = {"order": order, "graphs": graphs, "sha256": digest}
    missing = sorted(set(wanted) - set(selected))
    if missing:
        raise CheckViolation(
            f"{source}: no validated record for {missing}"
        )
    return selected


def select_r35_windows(path):
    """Complete R(3,5) edge windows, order zero synthesized as (0, 0)."""
    source = path.name
    raw = _load_json(path, source)
    results = raw.get("results")
    if not isinstance(results, list):
        raise CheckViolation(f"{source}: no results list")
    windows = {0: (0, 0)}       # the empty graph: no vertices, no edges
    for record in results:
        if not isinstance(record, dict):
            raise CheckViolation(f"{source}: malformed result record")
        if (record.get("s"), record.get("t")) != (3, 5):
            continue
        if not (record.get("ok") and record.get("count_matches")):
            continue            # unvalidated records cannot define a window
        order, _ = _record_order_and_count(record, "R(3,5) record", source)
        low = _exact_int(record.get("edge_min"), f"R(3,5,{order}) edge_min")
        high = _exact_int(record.get("edge_max"), f"R(3,5,{order}) edge_max")
        if low > high:
            raise CheckViolation(
                f"r35_edge_windows: order {order} window is empty "
                f"({low} > {high})"
            )
        known = windows.setdefault(order, (low, high))
        if known != (low, high):
            raise CheckViolation(
                f"r35_edge_windows: conflicting windows for order {order}: "
                f"{known} and {(low, high)}"
            )
    missing = [
        order for order in range(1, R35_MAX_ORDER + 1) if order not in windows
    ]
    if missing:
        raise CheckViolation(
            f"r35_edge_windows: no validated R(3,5) window for orders "
            f"{missing}"
        )
    return windows


def _extreme_name_class(name, source):
    """The (order, edges) class encoded by an `r45<order>.<edges>.g6` name."""
    if not name.startswith("r45") or not name.endswith(".g6"):
        raise CheckViolation(
            f"{source}: {name} is not an r45 edge-class catalog"
        )
    fields = name[3:-3].split(".")
    if len(fields) != 2 or not all(
        field.isascii() and field.isdigit() for field in fields
    ):
        raise CheckViolation(
            f"{source}: {name} does not encode order and edges"
        )
    return int(fields[0]), int(fields[1])


def select_extreme_records(path):
    """Every validated fixed-edge R(4,5) class at orders 17..23.

    Exactly 36 classes, every order of the band covered.  This report
    must carry no `count_matches` or `sha256` field: absence is
    required, not tolerated; filename, metadata order, and metadata
    edge class must agree instead.
    """
    source = path.name
    raw = _load_json(path, source)
    results = raw.get("results")
    if not isinstance(results, list):
        raise CheckViolation(f"{source}: no results list")
    selected = {}
    for record in results:
        if not isinstance(record, dict):
            raise CheckViolation(f"{source}: malformed result record")
        name, parent = _record_filename(record, source)
        order, graphs = _record_order_and_count(record, name, source)
        if order not in R45_EXTREME_ORDERS:
            continue            # the archive holds orders 4..16 as well
        if parent != R45_EXTREME_DIRNAME:
            raise CheckViolation(
                f"{source}: target-order record {name} is outside "
                f"{R45_EXTREME_DIRNAME}/"
            )
        if (record.get("s"), record.get("t")) != (4, 5):
            raise CheckViolation(f"{source}: {name} is not an R(4,5) record")
        _require_validated(record, ("ok", "all_ramsey"), name, source)
        for field in ("count_matches", "sha256"):
            if field in record:
                raise CheckViolation(
                    f"{source}: {name} carries a forbidden {field} field"
                )
        edge_min = _exact_int(record.get("edge_min"), f"{name} edge_min")
        edge_max = _exact_int(record.get("edge_max"), f"{name} edge_max")
        if edge_min != edge_max:
            raise CheckViolation(
                f"{source}: {name} spans edges {edge_min}..{edge_max}, not "
                "one fixed class"
            )
        if _extreme_name_class(name, source) != (order, edge_min):
            raise CheckViolation(
                f"{source}: {name} disagrees with its metadata "
                f"(order {order}, edges {edge_min})"
            )
        if name in selected:
            raise CheckViolation(
                f"duplicate: {source}: duplicate record for {name}"
            )
        selected[name] = {"order": order, "edges": edge_min, "graphs": graphs}
    if len(selected) != R45_EXTREME_RECORDS:
        raise CheckViolation(
            f"extreme: {len(selected)} validated order-"
            f"{R45_EXTREME_ORDERS[0]}..{R45_EXTREME_ORDERS[-1]} records in "
            f"{source}, expected {R45_EXTREME_RECORDS}"
        )
    covered = {info["order"] for info in selected.values()}
    missing = [order for order in R45_EXTREME_ORDERS if order not in covered]
    if missing:
        raise CheckViolation(
            f"extreme: no validated edge class for orders {missing}"
        )
    return selected


def load_edge_bounds(path):
    """Published R(4,5) edge extrema per order, from the owning table."""
    raw = _load_json(path, TABLE_FILE)
    table = raw.get("extremal_edge_tables_R45")
    if not isinstance(table, dict):
        raise CheckViolation(
            f"{TABLE_FILE}: no extremal_edge_tables_R45 object"
        )
    bounds = {}
    for order, values in table.items():
        if not isinstance(values, dict):
            raise CheckViolation(
                f"{TABLE_FILE}: order {order!r} is not an object"
            )
        low = _exact_int(values.get("min"), f"order {order} min")
        high = _exact_int(values.get("max"), f"order {order} max")
        if low > high:
            raise CheckViolation(
                f"{TABLE_FILE}: order {order} window is empty ({low} > {high})"
            )
        bounds[int(order)] = (low, high)
    missing = [order for order in STATE_DEGREES if order not in bounds]
    dual_missing = [
        order for order in (N_TARGET - 1 - d for d in STATE_DEGREES)
        if order not in bounds
    ]
    missing = sorted(set(missing) | set(dual_missing))
    if missing:
        raise CheckViolation(
            f"{TABLE_FILE}: no edge extrema for cone orders {missing}"
        )
    return bounds


# ------------------------------------------------- envelope and cone math ---

def degree_bounds(order):
    """Ramsey degree window of an R(4,5,order) graph."""
    return max(0, order - R44), min(R35_MAX_ORDER, order - 1)


def _iter_degree_histograms(low, high, vertices, degree_sum):
    """Every integer histogram over degrees low..high with fixed order/sum.

    Counts are indexed from `low`.  Graphicality is deliberately not
    imposed: the family is an outer relaxation, so skipping it only
    widens the envelope.  Depth-first with remaining-sum bounds.
    """
    if high < low or vertices < 0 or degree_sum < 0:
        return
    degrees = tuple(range(low, high + 1))
    top = degrees[-1]
    counts = [0] * len(degrees)

    def descend(position, verts_left, sum_left):
        degree = degrees[position]
        if position == len(degrees) - 1:
            if degree * verts_left == sum_left:
                counts[position] = verts_left
                yield tuple(counts)
                counts[position] = 0
            return
        floor_degree = degrees[position + 1]
        for count in range(verts_left + 1):
            sum_rest = sum_left - count * degree
            if sum_rest < 0:
                break           # larger counts only owe more
            verts_rest = verts_left - count
            if sum_rest > verts_rest * top:
                break           # the cap shrinks faster than the sum
            if sum_rest < verts_rest * floor_degree:
                continue        # too few vertices left to carry the sum
            counts[position] = count
            yield from descend(position + 1, verts_rest, sum_rest)
            counts[position] = 0

    yield from descend(0, vertices, degree_sum)


def q_outer_interval(order, edges, windows):
    """Outer interval for q = sum_x binom(d_x, 2) - t over an edge class.

    Every integer degree histogram of the class is enumerated; for each,
    the R(3,5) windows bound 3t to [T3lo, T3hi], integrality of t refines
    that to [ceil(T3lo/3), floor(T3hi/3)], and the local q interval is
    [W - t_hi, W - t_lo].  The envelope is the union over histograms.
    """
    if order < 0 or edges < 0:
        raise InfeasibleEdgeClass(
            f"negative class (order={order}, edges={edges})"
        )
    low, high = degree_bounds(order)
    degrees = range(low, high + 1)
    missing = [degree for degree in degrees if degree not in windows]
    if missing:
        raise CheckViolation(
            f"r35_edge_windows: no R(3,5) window for degrees {missing} "
            f"(order={order}, edges={edges})"
        )
    wedges = [degree * (degree - 1) // 2 for degree in degrees]
    tri3_lo = [windows[degree][0] for degree in degrees]
    tri3_hi = [windows[degree][1] for degree in degrees]
    q_lo = None
    q_hi = None
    standing = 0
    for histogram in _iter_degree_histograms(low, high, order, 2 * edges):
        standing += 1
        wedge_sum = 0
        lo_sum = 0
        hi_sum = 0
        for count, wedge, w_lo, w_hi in zip(
            histogram, wedges, tri3_lo, tri3_hi
        ):
            if count:
                wedge_sum += count * wedge
                lo_sum += count * w_lo
                hi_sum += count * w_hi
        t_lo = -(-lo_sum // 3)  # ceil towards the sound interior
        t_hi = hi_sum // 3
        if t_lo > t_hi:
            continue
        local_lo = wedge_sum - t_hi
        local_hi = wedge_sum - t_lo
        q_lo = local_lo if q_lo is None else min(q_lo, local_lo)
        q_hi = local_hi if q_hi is None else max(q_hi, local_hi)
    if q_lo is None:
        if standing:
            raise InfeasibleEdgeClass(
                f"no histogram meets the R(3,5) triangle window for "
                f"order={order}, edges={edges}"
            )
        raise InfeasibleEdgeClass(
            f"no degree histogram for order={order}, edges={edges}"
        )
    return q_lo, q_hi


def _resolve_interval(order, edges, windows, catalog_windows, cache):
    """Exact catalog window when the class is owned, else the relaxation.

    `None` means the class is infeasible: no integer degree histogram
    survives, so no state may use that side.
    """
    key = (order, edges)
    if key not in cache:
        if key in catalog_windows:
            cache[key] = catalog_windows[key]
        else:
            try:
                cache[key] = q_outer_interval(order, edges, windows)
            except InfeasibleEdgeClass:
                cache[key] = None
    return cache[key]


def budget2(d, edge_bounds):
    """Twice the largest m=2 excess deficit a degree-d vertex can carry."""
    m = N_TARGET - 1 - d
    h2_min = (
        2 * (m * (m - 1) // 2 - edge_bounds[m][1] - edge_bounds[d][1])
        - d * (N_TARGET - 2 * d)
    )
    return -h2_min


def _cone_state(d, a, b, ex, ey, qx, qy, budget):
    """One cone cell; the sole g formula with the frozen sign pairing."""
    m = N_TARGET - 1 - d
    dual_base = 3 * (m * (m - 1) * (m - 2) // 6 - ey * (m - 2))
    rhs_base = (N_TARGET + 3 - 3 * d) * ex
    return {
        "d": d,
        "a": a,
        "b": b,
        "deficiency": a + b,
        "excess_balance": 2 * (a + b) - budget,
        "g_lo": dual_base + 3 * qy[0] - rhs_base - 3 * qx[1],
        "g_hi": dual_base + 3 * qy[1] - rhs_base - 3 * qx[0],
    }


def rebuild_states(windows, edge_bounds, catalog_windows, class_source):
    """Every state of the n=45 cone, sorted by (d, a, b).

    A pair is omitted only when one side's edge class has no integer
    degree histogram at all; a missing catalog is never infeasibility,
    because such a class keeps its outer interval.  Per-side provenance
    strings are rebuilt from the observed class ownership.
    """
    cache = {}
    states = []
    for d in STATE_DEGREES:
        m = N_TARGET - 1 - d
        d_low, d_high = edge_bounds[d]
        m_low, m_high = edge_bounds[m]
        budget = budget2(d, edge_bounds)
        for a in range(d_high - d_low + 1):
            ex = d_high - a
            qx = _resolve_interval(
                d, ex, windows, catalog_windows, cache
            )
            if qx is None:
                continue
            for b in range(m_high - m_low + 1):
                ey = m_high - b
                qy = _resolve_interval(
                    m, ey, windows, catalog_windows, cache
                )
                if qy is None:
                    continue
                state = _cone_state(d, a, b, ex, ey, qx, qy, budget)
                state["x_interval_source"] = class_source.get(
                    (d, ex), OUTER_SOURCE
                )
                state["y_interval_source"] = class_source.get(
                    (m, ey), OUTER_SOURCE
                )
                states.append(state)
    states.sort(key=lambda state: (state["d"], state["a"], state["b"]))
    return states


# ---------------------------------------------------------- catalog rescan ---

def _stream_catalog(path, expected_order, expected_edges, expected_count):
    """One pass over a graph6 catalog: q frequencies per (order, edges).

    Every non-empty line is parsed, its order checked, its edge count
    checked when the class is fixed, and q accumulated.  No graph is
    retained; memory stays proportional to the tiny frequency maps.
    """
    label = path.name
    counts = {}
    graphs = 0
    try:
        stream = path.open("r", encoding="ascii")
    except OSError as exc:
        raise CheckViolation(
            f"catalog: {label}: unreadable catalog ({exc})"
        ) from None
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
                    raise CheckViolation(
                        f"catalog: {label} line {number}: bad graph6 ({exc})"
                    ) from None
                if order != expected_order:
                    raise CheckViolation(
                        f"catalog: {label} line {number}: order {order}, "
                        f"expected {expected_order}"
                    )
                wedges = 0
                degree_sum = 0
                for row in adj:
                    degree = popcount(row)
                    degree_sum += degree
                    wedges += degree * (degree - 1) // 2
                edges = degree_sum // 2
                if expected_edges is not None and edges != expected_edges:
                    raise CheckViolation(
                        f"catalog: {label} line {number}: {edges} edges, "
                        f"expected {expected_edges}"
                    )
                q = wedges - triangle_count(adj)
                frequency = counts.get((order, edges))
                if frequency is None:
                    frequency = counts[(order, edges)] = {}
                frequency[q] = frequency.get(q, 0) + 1
        except (OSError, UnicodeError) as exc:
            raise CheckViolation(
                f"catalog: {label}: unreadable catalog ({exc})"
            ) from None
    if not graphs:
        raise CheckViolation(f"catalog: {label}: empty catalog")
    if graphs != expected_count:
        raise CheckViolation(
            f"count: {label} holds {graphs} graphs, the validation report "
            f"records {expected_count}"
        )
    return graphs, counts


def sweep_catalogs(records):
    """Stream every catalog once: class histograms, line counts, owners.

    `records` is an iterable of `(path, expected_order,
    expected_edges_or_none, expected_graph_count, relative_path)`.  A
    class produced by two different files is an ownership violation.
    """
    histograms = {}
    owners = {}
    graph_counts = {}
    for path, order, edges, count, relative in records:
        location = path.as_posix()
        if location in graph_counts:
            raise CheckViolation(
                f"catalog: {path.name} is listed twice"
            )
        graphs, produced = _stream_catalog(path, order, edges, count)
        graph_counts[location] = graphs
        for key, frequency in produced.items():
            owner = owners.get(key)
            if owner is not None:
                raise CheckViolation(
                    f"owner: class (order={key[0]}, edges={key[1]}) is "
                    f"produced by both {owner} and {relative}"
                )
            owners[key] = relative
            histograms[key] = {
                q: frequency[q] for q in sorted(frequency)
            }
    return histograms, owners, graph_counts


def count_graph6_lines(path):
    """Non-empty graph6 line count of one small catalog."""
    count = 0
    try:
        with path.open("r", encoding="ascii") as stream:
            for line in stream:
                if line.strip():
                    count += 1
    except (OSError, UnicodeError) as exc:
        raise CheckViolation(
            f"count: {path.name}: unreadable catalog ({exc})"
        ) from None
    if not count:
        raise CheckViolation(f"count: {path.name}: empty catalog")
    return count


# ------------------------------------------------------------------ replay ---

def _check_replay(document, root):
    """Re-derive the replay section from the published file's own bytes.

    Every published graph and its complement must satisfy the checker's
    residual kernel exactly, and the recorded counts must be exactly
    what the file yields: lines, lines, and twice lines.
    """
    replay = document["replay"]
    path = root / "r55" / DATA_DIRNAME / R55_REPLAY_FILE
    name = path.name
    lines = 0
    zeros = 0
    try:
        stream = path.open("r", encoding="ascii")
    except OSError as exc:
        raise CheckViolation(
            f"replay: {name}: unreadable replay file ({exc})"
        ) from None
    with stream:
        try:
            for number, line in enumerate(stream, start=1):
                text = line.strip()
                if not text:
                    continue
                lines += 1
                try:
                    order, adj = parse_graph6_line(text)
                except ValueError as exc:
                    raise CheckViolation(
                        f"replay: {name} line {number}: bad graph6 ({exc})"
                    ) from None
                if order != R55_REPLAY_ORDER:
                    raise CheckViolation(
                        f"replay: {name} line {number}: order {order}, "
                        f"expected {R55_REPLAY_ORDER}"
                    )
                for label, graph in (
                    ("published graph", adj),
                    ("complement", complement(adj)),
                ):
                    residual = kernel_residual_scaled(graph)
                    if residual != 0:
                        raise CheckViolation(
                            f"residual: {name} line {number}: {label} has "
                            f"residual {residual}, expected 0"
                        )
                    zeros += 1
        except (OSError, UnicodeError) as exc:
            raise CheckViolation(
                f"replay: {name}: unreadable replay file ({exc})"
            ) from None
    if replay["published_graphs"] != lines:
        raise CheckViolation(
            f"replay: published_graphs is {replay['published_graphs']}, "
            f"the file holds {lines}"
        )
    if replay["complements"] != lines:
        raise CheckViolation(
            f"replay: complements is {replay['complements']}, the file "
            f"holds {lines}"
        )
    if replay["residual_zero"] != 2 * lines:
        raise CheckViolation(
            f"replay residual_zero: {replay['residual_zero']} recorded, "
            f"the replay verifies {2 * lines}"
        )
    return lines, zeros


# -------------------------------------------------------- search evidence ---

_ROUTES = {
    # objective_id: (per-state value, frozen acceptance edge, inclusive)
    "total_deficiency": (
        lambda state: state["deficiency"], Fraction(315), True,
    ),
    "degree20_count": (
        lambda state: int(state["d"] == 20), Fraction(1), False,
    ),
    "deficiency_ge8_count": (
        lambda state: int(state["deficiency"] >= 8), Fraction(1), False,
    ),
}


def _route(objective_id):
    """The frozen route of an objective; an unregistered ID is rejected."""
    route = _ROUTES.get(objective_id)
    if route is None:
        raise CheckViolation(
            f"objective_id: {objective_id!r} is not a frozen route"
        )
    return route


def _clears_edge(value, limit, inclusive):
    """True when an exact value clears the frozen acceptance edge."""
    return value <= limit if inclusive else value < limit


def _verify_certificate(states, entry, objective_id):
    """Check one affine certificate exactly; return the exact 45*alpha.

    Requires canonical coefficients with gamma, delta >= 0 and the
    inequality alpha + beta*xb + gamma*g_lo - delta*g_hi - objective >= 0
    at every single state.  Nothing here trusts how the coefficients
    were found.
    """
    _require_keys(entry, _CERTIFICATE_KEYS, "certificate")
    if entry["objective_id"] != objective_id:
        raise CheckViolation(
            f"certificate: objective {entry['objective_id']!r} is not the "
            f"record objective {objective_id!r}"
        )
    alpha = _canonical_fraction(entry["alpha"], "certificate alpha")
    beta = _canonical_fraction(entry["beta"], "certificate beta")
    gamma = _canonical_fraction(entry["gamma"], "certificate gamma")
    delta = _canonical_fraction(entry["delta"], "certificate delta")
    if gamma < 0:
        raise CheckViolation(f"certificate gamma: {gamma} is negative")
    if delta < 0:
        raise CheckViolation(f"certificate delta: {delta} is negative")
    value = _route(objective_id)[0]
    for index, state in enumerate(states):
        slack = (
            alpha
            + beta * state["excess_balance"]
            + gamma * state["g_lo"]
            - delta * state["g_hi"]
            - value(state)
        )
        if slack < 0:
            raise CheckViolation(
                f"certificate: state {index} (d={state['d']}, a={state['a']},"
                f" b={state['b']}) violates by {-slack} under alpha {alpha},"
                f" beta {beta}, gamma {gamma}, delta {delta}"
            )
    return VERTICES * alpha


def _verify_witness(states, entry, objective_id):
    """Check one primal witness exactly; return its exact objective value.

    Strictly increasing in-range indices, canonical nonnegative weights,
    total weight 45, zero aggregate excess balance, and aggregate
    g_lo <= 0 <= g_hi: exactly the relations a real 45-vertex graph
    satisfies.
    """
    _require_keys(entry, _WITNESS_KEYS, "primal witness")
    if entry["objective_id"] != objective_id:
        raise CheckViolation(
            f"primal witness: objective {entry['objective_id']!r} is not "
            f"the record objective {objective_id!r}"
        )
    weights = entry["weights"]
    if not isinstance(weights, list) or not weights:
        raise CheckViolation("primal witness: carries no weights")
    objective = _route(objective_id)[0]
    total = balance = low = high = value = Fraction(0)
    previous = None
    for item in weights:
        _require_keys(item, _WEIGHT_KEYS, "witness weight")
        index = _exact_int(item["state_index"], "witness state_index")
        if not 0 <= index < len(states):
            raise CheckViolation(
                f"state_index: {index} is out of the state range"
            )
        if previous is not None and index <= previous:
            raise CheckViolation(
                f"state_index: {index} does not increase past {previous}"
            )
        previous = index
        weight = _canonical_fraction(
            item["weight"], f"weight of state {index}"
        )
        if weight < 0:
            raise CheckViolation(
                f"weight: {weight} of state {index} is negative"
            )
        state = states[index]
        total += weight
        balance += weight * state["excess_balance"]
        low += weight * state["g_lo"]
        high += weight * state["g_hi"]
        value += weight * objective(state)
    if total != VERTICES:
        raise CheckViolation(
            f"weight: witness weights sum to {total}, not {VERTICES}"
        )
    if balance != 0:
        raise CheckViolation(
            f"witness: aggregate excess balance is {balance}, not zero"
        )
    if low > 0:
        raise CheckViolation(
            f"witness: aggregate g_lo is {low}, which exceeds zero"
        )
    if high < 0:
        raise CheckViolation(
            f"witness: aggregate g_hi is {high}, which is below zero"
        )
    return value


def _stored_exact(text, expected, label):
    """A stored exact value exists exactly when its evidence does."""
    if expected is None:
        if text is not None:
            raise CheckViolation(f"{label}: {text!r} has no evidence")
        return
    if text is None:
        raise CheckViolation(
            f"{label}: is null but the evidence gives {expected}"
        )
    if _canonical_fraction(text, label) != expected:
        raise CheckViolation(
            f"{label}: {text!r} is not the verified {expected}"
        )


def _verify_route_record(states, record):
    """Verify one route record against the artifact's own stored states."""
    _require_keys(record, _SEARCH_RECORD_KEYS, "search record")
    objective_id = record["objective_id"]
    objective, limit, inclusive = _route(objective_id)
    bound = None
    if record["certificate"] is not None:
        bound = _verify_certificate(
            states, record["certificate"], objective_id
        )
    _stored_exact(record["exact_upper_bound"], bound, "exact_upper_bound")
    value = None
    if record["primal_witness"] is not None:
        value = _verify_witness(
            states, record["primal_witness"], objective_id
        )
    _stored_exact(record["exact_witness_value"], value, "exact_witness_value")
    status = record["route_status"]
    if status not in ROUTE_STATUSES:
        raise CheckViolation(
            f"route_status: {status!r} is not a route status"
        )
    accepted = bound is not None and _clears_edge(bound, limit, inclusive)
    rejected = value is not None and not _clears_edge(
        value, limit, inclusive
    )
    if accepted and rejected:
        raise CheckViolation(
            f"route_status: accepted bound {bound} contradicts the witness "
            f"value {value} of {objective_id}"
        )
    if accepted:
        derived = ACCEPTED
    elif rejected:
        derived = REJECTED
    else:
        derived = UNRESOLVED
    if status != derived:
        raise CheckViolation(
            f"route_status: {status!r} is not the derived {derived} of "
            f"{objective_id}"
        )
    return derived


def _check_searches(document):
    """Verify all four search records; return the three route statuses."""
    searches = document["searches"]
    expected = len(ROUTE_IDS) + 1
    if not isinstance(searches, list):
        raise CheckViolation("searches: not a list")
    if len(searches) != expected:
        raise CheckViolation(
            f"searches: {len(searches)} records, expected exactly {expected}"
        )
    states = document["states"]
    statuses = []
    for index, objective_id in enumerate(ROUTE_IDS):
        record = searches[index]
        _require_keys(record, _SEARCH_RECORD_KEYS, "search record")
        if record["objective_id"] != objective_id:
            raise CheckViolation(
                f"route order: searches[{index}] objective "
                f"{record['objective_id']!r} is not the frozen route order "
                f"{objective_id!r}"
            )
        statuses.append(_verify_route_record(states, record))
    if searches[-1] != UNAVAILABLE_RECORD:
        raise CheckViolation(
            f"unavailable: searches[-1] must be the frozen "
            f"{UNAVAILABLE_ROUTE} record carrying no evidence, got "
            f"{searches[-1]!r}"
        )
    return statuses


def _terminal_disposition(statuses):
    """The verdict implied by the three route statuses alone."""
    if any(status == ACCEPTED for status in statuses):
        return DISPOSITION_ACCEPTED
    if all(status == REJECTED for status in statuses):
        return DISPOSITION_NO_CUT
    return DISPOSITION_UNRESOLVED


# ---------------------------------------------------------------- schema ---

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
_SEARCH_RECORD_KEYS = (
    "certificate", "exact_upper_bound", "exact_witness_value", "objective_id",
    "primal_witness", "route_status",
)
_CERTIFICATE_KEYS = ("alpha", "beta", "delta", "gamma", "objective_id")
_WITNESS_KEYS = ("objective_id", "weights")
_WEIGHT_KEYS = ("state_index", "weight")


def _check_relative_path(record):
    """A canonical `r55/...` POSIX path; nothing absolute may appear."""
    label = record["relative_path"]
    if type(label) is not str or not label:
        raise CheckViolation("relative_path: input record without a path")
    pure = PurePosixPath(label)
    if pure.is_absolute() or label.startswith("/"):
        raise CheckViolation(f"relative_path: {label!r} is absolute")
    parts = pure.parts
    if not parts or parts[0] != "r55":
        raise CheckViolation(
            f"relative_path: {label!r} is not rooted at r55/"
        )
    if any(part in ("", ".", "..") for part in parts):
        raise CheckViolation(f"relative_path: {label!r} is not canonical")
    return label


def _check_schema(document):
    """Exact key sets and per-section self-consistency, no disk access."""
    if not isinstance(document, dict):
        raise CheckViolation(f"schema: document is not an object")
    _require_keys(document, _TOP_LEVEL_KEYS, "schema: analysis document")
    if _exact_int(document["schema_version"], "schema_version") != (
        SCHEMA_VERSION
    ):
        raise CheckViolation(
            f"schema_version: {document['schema_version']!r}"
        )
    if document["campaign_id"] != CAMPAIGN_ID:
        raise CheckViolation(f"campaign_id: {document['campaign_id']!r}")
    if type(document["disposition"]) is not str:
        raise CheckViolation(
            f"disposition: {document['disposition']!r} is not a string"
        )
    for text in _iter_strings(document):
        if text.startswith("/"):
            raise CheckViolation(
                f"schema: absolute path {text!r} in the analysis"
            )

    inputs = document["inputs"]
    if not isinstance(inputs, list) or not inputs:
        raise CheckViolation("inputs: no input records")
    paths = []
    for record in inputs:
        _require_keys(record, _INPUT_KEYS, "input record")
        if not _is_sha256(record["sha256"]):
            raise CheckViolation(
                f"sha256: {record.get('relative_path')!r} is malformed"
            )
        if _exact_int(record["bytes"], "bytes") <= 0:
            raise CheckViolation("bytes: empty input")
        count = record["graph_count"]
        if count is not None and _exact_int(count, "graph_count") <= 0:
            raise CheckViolation("graph_count: non-positive")
        paths.append(_check_relative_path(record))
    if paths != sorted(set(paths)):
        raise CheckViolation(
            "inputs: records are not sorted and unique by relative_path"
        )

    roots = document["trust_roots"]
    if not isinstance(roots, list) or not roots:
        raise CheckViolation("trust_roots: must be a non-empty list")
    for record in roots:
        _require_keys(record, _TRUST_KEYS, "trust root")
        if type(record["statement"]) is not str:
            raise CheckViolation(
                f"trust_roots: statement {record['statement']!r} is not "
                "a string"
            )
        if not record["statement"]:
            raise CheckViolation("trust_roots: root without a statement")
        if type(record["url"]) is not str:
            raise CheckViolation(
                f"trust_roots: url {record['url']!r} is not a string"
            )
        if not record["url"].startswith("https://"):
            raise CheckViolation(
                f"trust_roots: url {record['url']!r} is not https"
            )
    if roots != sorted(roots, key=lambda r: (r["statement"], r["url"])):
        raise CheckViolation("trust_roots: not sorted")

    replay = document["replay"]
    _require_keys(replay, _REPLAY_KEYS, "replay")
    for key in _REPLAY_KEYS:
        _exact_int(replay[key], f"replay {key}")

    r35_windows = document["r35_edge_windows"]
    if not isinstance(r35_windows, list):
        raise CheckViolation("r35_edge_windows: not a list")
    orders = []
    for record in r35_windows:
        _require_keys(record, _R35_WINDOW_KEYS, "R(3,5) window")
        order = _exact_int(record["order"], "R(3,5) window order")
        low = _exact_int(record["edge_min"], "R(3,5) window edge_min")
        high = _exact_int(record["edge_max"], "R(3,5) window edge_max")
        if low > high:
            raise CheckViolation(
                f"r35_edge_windows: order {order} window is empty"
            )
        orders.append(order)
    if orders != list(range(R35_MAX_ORDER + 1)):
        raise CheckViolation(
            f"r35_edge_windows: orders {orders} do not cover exactly "
            f"0..{R35_MAX_ORDER}"
        )

    histograms = document["catalog_q_histograms"]
    if not isinstance(histograms, list) or not histograms:
        raise CheckViolation("catalog_q_histograms: no classes")
    if len(histograms) <= R45_EXTREME_RECORDS:
        raise CheckViolation(
            f"catalog_q_histograms: {len(histograms)} classes, expected "
            f"more than the {R45_EXTREME_RECORDS} extreme classes"
        )
    observed = {}
    classes = []
    sources = set()
    for record in histograms:
        _require_keys(record, _HISTOGRAM_KEYS, "catalog histogram")
        order = _exact_int(record["order"], "catalog histogram order")
        edges = _exact_int(record["edges"], "catalog histogram edges")
        key = (order, edges)
        classes.append(key)
        source = record["source"]
        if not isinstance(source, str) or not source.endswith(".g6"):
            raise CheckViolation(
                f"source: class {key} has no graph6 source"
            )
        sources.add(source)
        entries = record["q_counts"]
        if not isinstance(entries, list) or not entries:
            raise CheckViolation(
                f"catalog_q_histograms: class {key} has an empty histogram"
            )
        values = []
        total = 0
        for entry in entries:
            _require_keys(entry, _Q_COUNT_KEYS, "q count")
            frequency = _exact_int(entry["count"], "q count count")
            if frequency <= 0:
                raise CheckViolation(
                    f"catalog_q_histograms: class {key} has a non-positive "
                    "frequency"
                )
            values.append(_exact_int(entry["q"], "q count q"))
            total += frequency
        if values != sorted(set(values)):
            raise CheckViolation(
                f"catalog_q_histograms: class {key} has unsorted or "
                "repeated q values"
            )
        graph_count = _exact_int(
            record["graph_count"], "catalog histogram graph_count"
        )
        if total != graph_count or graph_count <= 0:
            raise CheckViolation(
                f"graph_count: class {key} sums to {total}, not its "
                f"{graph_count} graphs"
            )
        observed[key] = (values[0], values[-1], source)
    if classes != sorted(set(classes)):
        raise CheckViolation(
            "catalog_q_histograms: classes are not sorted and unique"
        )

    exact_windows = document["catalog_windows"]
    if not isinstance(exact_windows, list):
        raise CheckViolation("catalog_windows: not a list")
    if len(exact_windows) != len(classes):
        raise CheckViolation(
            f"catalog_windows: {len(exact_windows)} windows for "
            f"{len(classes)} classes"
        )
    seen = []
    for record in exact_windows:
        _require_keys(record, _CATALOG_WINDOW_KEYS, "catalog window")
        order = _exact_int(record["order"], "catalog window order")
        edges = _exact_int(record["edges"], "catalog window edges")
        key = (order, edges)
        if key not in observed:
            raise CheckViolation(
                f"catalog_windows: window {key} has no histogram"
            )
        q_min = _exact_int(record["q_min"], "catalog window q_min")
        q_max = _exact_int(record["q_max"], "catalog window q_max")
        if (q_min, q_max, record["source"]) != observed[key]:
            raise CheckViolation(
                f"catalog_windows: window {key} disagrees with its own "
                "histogram"
            )
        seen.append(key)
    if seen != sorted(set(seen)):
        raise CheckViolation(
            "catalog_windows: classes are not sorted and unique"
        )

    states = document["states"]
    if not isinstance(states, list) or not states:
        raise CheckViolation("states: the analysis records no states")
    cells = []
    for record in states:
        _require_keys(record, _STATE_KEYS, "state")
        for field in ("a", "b", "d", "deficiency", "excess_balance",
                      "g_hi", "g_lo"):
            _exact_int(record[field], f"state {field}")
        if record["d"] not in STATE_DEGREES:
            raise CheckViolation(
                f"d: state degree {record['d']} outside "
                f"{STATE_DEGREES[0]}..{STATE_DEGREES[-1]}"
            )
        if record["a"] < 0 or record["b"] < 0:
            raise CheckViolation(
                f"state: negative shortfall at d={record['d']}"
            )
        if record["deficiency"] != record["a"] + record["b"]:
            raise CheckViolation(
                "deficiency: disagrees with a + b"
            )
        if record["g_lo"] > record["g_hi"]:
            raise CheckViolation("g_lo: state has an empty g interval")
        for field in ("x_interval_source", "y_interval_source"):
            value = record[field]
            if value != OUTER_SOURCE and value not in sources:
                raise CheckViolation(
                    f"{field}: {value!r} is not a catalog source"
                )
        cells.append((record["d"], record["a"], record["b"]))
    if cells != sorted(set(cells)):
        raise CheckViolation("states: not sorted and unique")

    searches = document["searches"]
    if not isinstance(searches, list) or len(searches) != (
        len(ROUTE_IDS) + 1
    ):
        raise CheckViolation(
            f"searches: expected exactly {len(ROUTE_IDS) + 1} records"
        )


# ---------------------------------------------------------------- hashes ---

def _check_hashes(document, root):
    """Every recorded input still hashes to its recorded SHA-256."""
    for record in document["inputs"]:
        label = record["relative_path"]
        path = root / label
        if not path.is_file():
            raise CheckViolation(
                f"inputs: {label} is missing under the root {root}"
            )
        digest = _sha256_of(path)
        if digest != record["sha256"]:
            raise CheckViolation(
                f"sha256: {label} records {record['sha256']}, the local "
                f"bytes hash to {digest}"
            )


def _check_sizes(document, root):
    """Every recorded byte count equals the file size on disk.

    Redundant with the digest in isolation, but it pins the recorded
    `bytes` field itself; it runs after the validation selection so a
    structurally broken report is named before any byte count is.
    """
    for record in document["inputs"]:
        label = record["relative_path"]
        try:
            size = (root / label).stat().st_size
        except OSError as exc:
            raise CheckViolation(
                f"bytes: {label}: unstatable input ({exc})"
            ) from None
        if size != record["bytes"]:
            raise CheckViolation(
                f"bytes: {label} records {record['bytes']}, the local file "
                f"holds {size}"
            )


# ---------------------------------------------------------- catalog stage ---

def _derived_inputs(main_records, extreme_records):
    """The required input set: rel path -> required graph_count or None."""
    required = {}
    for order, name in enumerate(R35_FILES, start=1):
        required[f"r55/{DATA_DIRNAME}/{name}"] = (
            main_records[name]["graphs"]
        )
    for name in sorted(extreme_records):
        required[
            f"r55/{DATA_DIRNAME}/{R45_EXTREME_DIRNAME}/{name}"
        ] = extreme_records[name]["graphs"]
    required[f"r55/{DATA_DIRNAME}/{R45_CENSUS_FILE}"] = (
        main_records[R45_CENSUS_FILE]["graphs"]
    )
    required[f"r55/{DATA_DIRNAME}/{R55_REPLAY_FILE}"] = (
        main_records[R55_REPLAY_FILE]["graphs"]
    )
    for name in (VALIDATION_FILE, VALIDATION_EXTREME_FILE, TABLE_FILE):
        required[f"r55/{DATA_DIRNAME}/{name}"] = None
    for name in RECORDED_SOURCES:
        required[f"r55/src/{name}"] = None
    return required


def _check_inputs_vs_derived(document, required):
    """The recorded inputs are exactly the derived required set."""
    recorded = {
        record["relative_path"]: record["graph_count"]
        for record in document["inputs"]
    }
    missing = sorted(set(required) - set(recorded))
    if missing:
        raise CheckViolation(f"inputs: no record for {missing}")
    extra = sorted(set(recorded) - set(required))
    if extra:
        raise CheckViolation(f"inputs: unexpected records {extra}")
    for label in sorted(required):
        if recorded[label] != required[label]:
            raise CheckViolation(
                f"graph_count: {label} records {recorded[label]}, the "
                f"validated files give {required[label]}"
            )


def _check_catalogs(document, root):
    """Re-select, re-count, and re-sweep every catalog the artifact claims.

    Returns the derived (windows, class_source, catalog_windows) for the
    state rebuild.
    """
    data = root / "r55" / DATA_DIRNAME
    main_records = select_main_records(data / VALIDATION_FILE)
    windows = select_r35_windows(data / VALIDATION_FILE)
    extreme_records = select_extreme_records(data / VALIDATION_EXTREME_FILE)

    _check_inputs_vs_derived(
        document, _derived_inputs(main_records, extreme_records)
    )

    _check_sizes(document, root)

    for order, name in enumerate(R35_FILES, start=1):
        graphs = count_graph6_lines(data / name)
        if graphs != main_records[name]["graphs"]:
            raise CheckViolation(
                f"count: {name} holds {graphs} non-empty lines, the "
                f"validation report records "
                f"{main_records[name]['graphs']}"
            )

    replay_lines = count_graph6_lines(data / R55_REPLAY_FILE)
    if replay_lines != main_records[R55_REPLAY_FILE]["graphs"]:
        raise CheckViolation(
            f"count: {R55_REPLAY_FILE} holds {replay_lines} non-empty "
            f"lines, the validation report records "
            f"{main_records[R55_REPLAY_FILE]['graphs']}"
        )

    recorded_shas = {
        record["relative_path"]: record["sha256"]
        for record in document["inputs"]
    }
    for name in (*R35_FILES, R45_CENSUS_FILE, R55_REPLAY_FILE):
        label = f"r55/{DATA_DIRNAME}/{name}"
        validated = main_records[name]["sha256"]
        if recorded_shas[label] != validated:
            raise CheckViolation(
                f"sha256: {label}: the artifact records "
                f"{recorded_shas[label]}, the validation report records "
                f"{validated}"
            )

    records = []
    for name in sorted(extreme_records):
        info = extreme_records[name]
        records.append((
            data / R45_EXTREME_DIRNAME / name,
            info["order"], info["edges"], info["graphs"],
            f"r55/{DATA_DIRNAME}/{R45_EXTREME_DIRNAME}/{name}",
        ))
    census_rel = f"r55/{DATA_DIRNAME}/{R45_CENSUS_FILE}"
    records.append((
        data / R45_CENSUS_FILE,
        main_records[R45_CENSUS_FILE]["order"], None,
        main_records[R45_CENSUS_FILE]["graphs"], census_rel,
    ))
    histograms, owners, _counts = sweep_catalogs(records)

    for key, owner in sorted(owners.items()):
        order, edges = key
        if owner == census_rel:
            if order != CENSUS_ORDER or not (
                CENSUS_EDGE_BAND[0] <= edges <= CENSUS_EDGE_BAND[1]
            ):
                raise CheckViolation(
                    f"owner: the census owns (order={order}, edges={edges}), "
                    f"outside the order-{CENSUS_ORDER} band "
                    f"{CENSUS_EDGE_BAND}"
                )
        elif order == CENSUS_ORDER:
            raise CheckViolation(
                f"owner: order-{CENSUS_ORDER} class {key} is owned by "
                f"{owner}, not the census"
            )

    doc_histograms = document["catalog_q_histograms"]
    doc_classes = [(r["order"], r["edges"]) for r in doc_histograms]
    observed_classes = sorted(histograms)
    if doc_classes != observed_classes:
        raise CheckViolation(
            f"catalog_q_histograms: recorded classes {doc_classes} are not "
            f"the observed {observed_classes}"
        )
    for record in doc_histograms:
        key = (record["order"], record["edges"])
        frequency = histograms[key]
        expected_counts = [
            {"q": q, "count": frequency[q]} for q in sorted(frequency)
        ]
        if record["q_counts"] != expected_counts:
            raise CheckViolation(
                f"q_counts: class {key} does not match the reswept "
                "frequencies"
            )
        if record["source"] != owners[key]:
            raise CheckViolation(
                f"source: class (order={key[0]}, edges={key[1]}) records "
                f"{record['source']!r}, the observed owner is "
                f"{owners[key]!r}"
            )

    for record in document["catalog_windows"]:
        key = (record["order"], record["edges"])
        frequency = histograms[key]
        if record["q_min"] != min(frequency):
            raise CheckViolation(
                f"q_min: class {key} records {record['q_min']}, observed "
                f"{min(frequency)}"
            )
        if record["q_max"] != max(frequency):
            raise CheckViolation(
                f"q_max: class {key} records {record['q_max']}, observed "
                f"{max(frequency)}"
            )
        if record["source"] != owners[key]:
            raise CheckViolation(
                f"source: class {key} window source is not the observed "
                f"owner {owners[key]!r}"
            )

    for key in observed_classes:
        try:
            low, high = q_outer_interval(key[0], key[1], windows)
        except InfeasibleEdgeClass as exc:
            raise CheckViolation(
                f"envelope: class (order={key[0]}, edges={key[1]}) has no "
                f"surviving degree histogram ({exc})"
            ) from None
        outside = [q for q in histograms[key] if not low <= q <= high]
        if outside:
            raise CheckViolation(
                f"envelope: class (order={key[0]}, edges={key[1]}) observes "
                f"q values {outside} outside its outer interval "
                f"[{low}, {high}]"
            )

    catalog_windows = {
        key: (min(histograms[key]), max(histograms[key]))
        for key in observed_classes
    }
    class_source = dict(owners)
    total = sum(sum(frequency.values()) for frequency in histograms.values())
    return windows, class_source, catalog_windows, total


# ------------------------------------------------------------ state stage ---

def _check_states(document, root, windows, class_source, catalog_windows):
    """Rebuild the whole cone and require every recorded field to match."""
    edge_bounds = load_edge_bounds(
        root / "r55" / DATA_DIRNAME / TABLE_FILE
    )

    derived_windows = [
        {
            "order": order,
            "edge_min": windows[order][0],
            "edge_max": windows[order][1],
        }
        for order in sorted(windows)
    ]
    if document["r35_edge_windows"] != derived_windows:
        raise CheckViolation(
            "r35_edge_windows: the recorded windows are not the validated "
            f"ones ({derived_windows})"
        )

    rebuilt = rebuild_states(windows, edge_bounds, catalog_windows,
                             class_source)
    states = document["states"]
    if len(states) != len(rebuilt):
        raise CheckViolation(
            f"states: {len(states)} recorded, the rebuild gives "
            f"{len(rebuilt)}"
        )
    for index, (record, mine) in enumerate(zip(states, rebuilt)):
        for field in ("d", "a", "b", "deficiency", "excess_balance",
                      "g_lo", "g_hi", "x_interval_source",
                      "y_interval_source"):
            if record[field] != mine[field]:
                raise CheckViolation(
                    f"{field}: states[{index}] (d={record['d']}, "
                    f"a={record['a']}, b={record['b']}) records "
                    f"{record[field]!r}, the rebuild gives {mine[field]!r}"
                )


# ---------------------------------------------------------------- verdict ---

def verify_analysis(document, root=None, echo=None):
    """Re-verify one analysis document end to end; return the disposition.

    Every stage fails fast with a `CheckViolation` naming the first
    broken field.  `root` substitutes the math repository root that all
    recorded `r55/...` paths resolve against; it defaults to the
    checker's own repository root.  `echo`, when given, receives one
    deterministic progress line per completed stage.
    """
    progress = echo if callable(echo) else (lambda line: None)
    if root is None:
        root = _DEFAULT_ROOT
    root = Path(root).resolve()

    _check_schema(document)
    _check_hashes(document, root)
    progress(f"M3 HASHES: {len(document['inputs'])} inputs byte-identical")

    _check_replay(document, root)
    progress(
        f"M3 REPLAY: {document['replay']['residual_zero']} of "
        f"{2 * document['replay']['published_graphs']} residuals zero"
    )

    statuses = _check_searches(document)
    progress(
        f"M3 SEARCHES: {len(ROUTE_IDS)} frozen routes verified against "
        f"{len(document['states'])} stored states"
    )

    windows, class_source, catalog_windows, total = _check_catalogs(
        document, root
    )
    progress(
        f"M3 CATALOGS: {total} graphs reswept into {len(class_source)} "
        "classes"
    )

    _check_states(document, root, windows, class_source, catalog_windows)
    progress(f"M3 STATES: {len(document['states'])} states rebuilt exactly")

    disposition = _terminal_disposition(statuses)
    if document["disposition"] != disposition:
        raise CheckViolation(
            f"disposition: {document['disposition']!r} is not the terminal "
            f"{disposition!r}"
        )
    _check_hashes(document, root)
    return disposition


def main(argv=None):
    """Verify one analysis artifact; print the verdict, return the code."""
    parser = argparse.ArgumentParser(
        description=(
            "Independently re-verify an m=3 higher-identity analysis "
            "artifact against its recorded inputs."
        )
    )
    parser.add_argument(
        "--root",
        default=None,
        help=(
            "math repository root every recorded r55/... path resolves "
            "against (default: the checker's own repository root)"
        ),
    )
    parser.add_argument(
        "analysis",
        help="path to higher_identity_m3.json",
    )
    args = parser.parse_args(argv)
    root = Path(args.root) if args.root is not None else _DEFAULT_ROOT
    try:
        target = Path(args.analysis)
        try:
            document = json.loads(target.read_text(encoding="ascii"))
        except (OSError, UnicodeError, json.JSONDecodeError) as exc:
            raise CheckViolation(
                f"analysis: {target}: unreadable JSON ({exc})"
            ) from None
        disposition = verify_analysis(document, root, echo=_echo)
    except CheckViolation as exc:
        print(f"M3 EVIDENCE FAILED: {exc}", file=sys.stderr)
        return 1
    except (OSError, ValueError) as exc:
        print(f"M3 EVIDENCE FAILED: analysis: {exc}", file=sys.stderr)
        return 1
    print(f"M3 EVIDENCE VERIFIED: {disposition}")
    return 0


def _echo(line):
    print(line, flush=True)


if __name__ == "__main__":
    sys.exit(main())
