"""Focused tests for the sound m=3 motif envelopes and the n=45 state cone.

Every oracle here is independent of `m3_deficiency_cone`: degree windows are
restated from the Ramsey numbers R(3,5)=14 and R(4,4)=18, small histogram
families are brute forced, exact motif counts come from the Task 1 module,
and the m=2 budget is cross-checked against `structural_constraints.h_bounds`
(the code that verified the excess identity on all 656 published graphs).
"""

import dataclasses
import itertools
import json
import pathlib
import sys
import tempfile
import unittest
from fractions import Fraction
from unittest import mock

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import structural_constraints
from check_ramsey import encode_graph6, parse_graph6_line
from m3_deficiency_cone import (  # expected to fail before Task 2 implementation
    InfeasibleEdgeClass, State, budget2, build_states, catalog_q_histograms,
    degree_histograms, load_edge_bounds, m3_state_interval, q_interval,
    q_outer_interval, r35_edge_windows,
)
from subgraph_identities import complement, m3_residual_scaled, motif3

DATA = ROOT / "data"
EXTREME = DATA / "r45extreme"
STATE_DEGREES = (20, 21, 22, 23, 24)
_SAMPLES = {}


def degree_window(order):
    """Ramsey degree window of an R(4,5,order) graph, restated from scratch.

    N(x) is K3-free and I5-free, hence an R(3,5) graph, so d_x <= 13; the
    graph minus N[x] is K4-free and I4-free, hence an R(4,4) graph, so
    order-1-d_x <= 17.
    """
    return max(0, order - 18), min(13, order - 1)


def read_graphs(path, edges=None, limit=None):
    """Graphs of a graph6 file, restricted to a given edge count if asked."""
    graphs = []
    with path.open(encoding="ascii") as handle:
        for line in handle:
            text = line.strip()
            if not text:
                continue
            if edges is not None and sum(
                (ord(char) - 63).bit_count() for char in text[1:]
            ) != edges:
                continue
            _, adj = parse_graph6_line(text)
            if edges is not None:
                assert sum(row.bit_count() for row in adj) == 2 * edges
            graphs.append(adj)
            if limit is not None and len(graphs) == limit:
                break
    return graphs


def catalog_sample(path, edges, limit):
    """Cached `read_graphs`; the catalogs are immutable inputs."""
    key = (path, edges, limit)
    if key not in _SAMPLES:
        _SAMPLES[key] = read_graphs(path, edges=edges, limit=limit)
    return _SAMPLES[key]


def q_exact(adj):
    """q(H) = sum_x binom(d_x, 2) - t(H), computed from Task 1 counts."""
    f = motif3(adj)
    return sum(
        row.bit_count() * (row.bit_count() - 1) // 2 for row in adj
    ) - f.triangles


def brute_force_histograms(order, edges):
    """All integer degree histograms, by exhaustive product enumeration."""
    low, high = degree_window(order)
    found = []
    for hist in itertools.product(range(order + 1), repeat=high - low + 1):
        if sum(hist) != order:
            continue
        if sum((low + i) * c for i, c in enumerate(hist)) == 2 * edges:
            found.append(hist)
    return found


def partitions(total, largest):
    """Every partition of `total` into descending parts of size 1..largest."""
    if total == 0:
        yield ()
        return
    for part in range(min(total, largest), 0, -1):
        for rest in partitions(total - part, part):
            yield (part,) + rest


def clamped_histograms(order, edges):
    """All integer degree histograms, via the Ramsey-clamp bijection.

    Writing d_x = low + y_x turns the degree window and the handshake
    equation into "partitions of 2*edges - order*low into at most `order`
    parts, each at most high - low". Enumerating those partitions is bounded
    by p(excess degree sum), where the product oracle above would need
    (order + 1) ** (high - low + 1) tuples. Validated against that oracle at
    orders 5 and 6 before it is trusted at order 20.
    """
    low, high = degree_window(order)
    found = []
    for parts in partitions(2 * edges - order * low, high - low):
        if len(parts) > order:
            continue
        counts = [0] * (high - low + 1)
        counts[0] = order - len(parts)
        for part in parts:
            counts[part] += 1
        found.append(tuple(counts))
    return found


def unit_windows():
    """Trivial R(3,5) windows: every k-vertex graph has 0..binom(k,2) edges."""
    return {k: (0, k * (k - 1) // 2) for k in range(14)}


class M3ConeTests(unittest.TestCase):
    def test_degree_histograms_obey_order_and_handshake(self):
        for hist in degree_histograms(5, 5):
            lo = max(0, 5 - 18)
            self.assertEqual(sum(hist), 5)
            self.assertEqual(sum((lo + i) * c for i, c in enumerate(hist)), 10)

    def test_degree_histograms_are_exhaustive_and_lexicographic(self):
        for order, edges in ((5, 5), (6, 7)):
            produced = list(degree_histograms(order, edges))
            self.assertEqual(produced, sorted(brute_force_histograms(order, edges)))
            self.assertEqual(
                sorted(clamped_histograms(order, edges)),
                sorted(brute_force_histograms(order, edges)),
            )
            self.assertEqual(produced, sorted(produced))

    def test_degree_histograms_respect_both_ramsey_clamps(self):
        # order 20 clamps below (20-18=2) and above (R(3,5)-1=13).
        low, high = degree_window(20)
        self.assertEqual((low, high), (2, 13))
        produced = []
        for hist in degree_histograms(20, 30):
            self.assertEqual(len(hist), high - low + 1)
            self.assertEqual(sum(hist), 20)
            self.assertEqual(
                sum((low + i) * c for i, c in enumerate(hist)), 60
            )
            produced.append(hist)
        # The product oracle would need 21 ** 12 tuples at order 20, so the
        # bounded clamp bijection pins the family here; orders 5 and 6 keep
        # the generator honest against the product oracle itself.
        self.assertEqual(produced, sorted(clamped_histograms(20, 30)))
        # Both clamps are genuinely exercised: the lower clamp 20-18=2 and
        # the upper clamp R(3,5)-1=13 both carry vertices, nothing outside.
        self.assertEqual(
            {low + i for hist in produced for i, c in enumerate(hist) if c},
            set(range(low, high + 1)),
        )

    def test_outer_interval_contains_exact_c5_value(self):
        windows = unit_windows()
        lo, hi = q_outer_interval(5, 5, windows)
        exact_q = 5  # five degree-2 wedges, zero triangles
        self.assertLessEqual(lo, exact_q)
        self.assertLessEqual(exact_q, hi)

    def test_exact_catalog_window_overrides_outer_relaxation(self):
        windows = unit_windows()
        catalog_windows = {(5, 5): (7, 7)}
        self.assertEqual(q_interval(5, 5, windows, catalog_windows), (7, 7))

    def test_outer_interval_covers_every_r35_catalog_graph(self):
        windows = r35_edge_windows(DATA / "VALIDATION.json")
        for order, edges in ((11, 22), (12, 24), (13, 26)):
            lo, hi = q_outer_interval(order, edges, windows)
            graphs = catalog_sample(DATA / f"r35_{order}.g6", edges, None)
            self.assertTrue(graphs)
            for adj in graphs:
                self.assertLessEqual(lo, q_exact(adj))
                self.assertLessEqual(q_exact(adj), hi)

    def test_state_interval_contains_exact_catalog_pairs(self):
        windows = r35_edge_windows(DATA / "VALIDATION.json")
        edge_bounds = load_edge_bounds(DATA / "structural_tables.json")
        max_edges = {order: bounds[1] for order, bounds in edge_bounds.items()}
        pairs = (
            (20, 100, EXTREME / "r4520.100.g6", 24, 132, DATA / "r45_24.g6"),
            (21, 107, EXTREME / "r4521.107.g6",
             23, 122, EXTREME / "r4523.122.g6"),
            (22, 114, EXTREME / "r4522.114.g6",
             22, 114, EXTREME / "r4522.114.g6"),
        )
        for d, ex, x_path, m, ey, y_path in pairs:
            self.assertEqual(m, 44 - d)
            state = m3_state_interval(
                d, max_edges[d] - ex, max_edges[m] - ey, windows, edge_bounds
            )
            xs = catalog_sample(x_path, ex, 20)
            ys = catalog_sample(y_path, ey, 20)
            self.assertTrue(xs)
            self.assertTrue(ys)
            for adj in xs:
                self.assertEqual(len(adj), d)
            for adj in ys:
                self.assertEqual(len(adj), m)
            for qx, qy in itertools.product(
                [q_exact(adj) for adj in xs], [q_exact(adj) for adj in ys]
            ):
                exact = (
                    3 * (m * (m - 1) * (m - 2) // 6 - ey * (m - 2) + qy)
                    - (48 - 3 * d) * ex - 3 * qx
                )
                self.assertLessEqual(state.g_lo, exact)
                self.assertLessEqual(exact, state.g_hi)

    def test_all_656_known_r55_graphs_have_zero_m3_residual(self):
        published = read_graphs(DATA / "r55_42some.g6")
        self.assertEqual(len(published), 328)
        tested = 0
        for adj in published:
            self.assertEqual(m3_residual_scaled(adj), 0)
            self.assertEqual(m3_residual_scaled(complement(adj)), 0)
            tested += 2
        self.assertEqual(tested, 656)

    def test_r35_edge_windows_match_validation_and_synthesize_order_zero(self):
        windows = r35_edge_windows(DATA / "VALIDATION.json")
        self.assertEqual(sorted(windows), list(range(14)))
        self.assertEqual(windows[0], (0, 0))
        self.assertEqual(windows[13], (26, 26))
        records = json.loads((DATA / "VALIDATION.json").read_text())["results"]
        checked = 0
        for record in records:
            if (record["s"], record["t"]) != (3, 5):
                continue
            order = int(next(iter(record["n_values"])))
            self.assertEqual(
                windows[order], (record["edge_min"], record["edge_max"])
            )
            checked += 1
        self.assertEqual(checked, 13)

    def test_load_edge_bounds_matches_structural_tables(self):
        edge_bounds = load_edge_bounds(DATA / "structural_tables.json")
        table = json.loads(
            (DATA / "structural_tables.json").read_text()
        )["extremal_edge_tables_R45"]
        self.assertEqual(sorted(edge_bounds), sorted(int(k) for k in table))
        for order, (low, high) in edge_bounds.items():
            self.assertEqual((low, high),
                             (table[str(order)]["min"], table[str(order)]["max"]))
            self.assertIs(type(low), int)
            self.assertIs(type(high), int)
            self.assertLessEqual(low, high)

    def test_budget2_matches_independent_m2_h_bound(self):
        edge_bounds = load_edge_bounds(DATA / "structural_tables.json")
        minima = {order: bounds[0] for order, bounds in edge_bounds.items()}
        maxima = {order: bounds[1] for order, bounds in edge_bounds.items()}
        with mock.patch.dict(structural_constraints.e45, minima), \
                mock.patch.dict(structural_constraints.E45, maxima):
            for d in STATE_DEGREES:
                h2_min = 2 * structural_constraints.h_bounds(45, d)[0]
                self.assertEqual(h2_min.denominator, 1)
                self.assertEqual(
                    Fraction(budget2(d, edge_bounds)), -h2_min
                )

    def test_budget2_reproduces_published_n45_deficiency_budget(self):
        raw = json.loads((DATA / "structural_tables.json").read_text())
        edge_bounds = load_edge_bounds(DATA / "structural_tables.json")
        worst = max(budget2(d, edge_bounds) for d in STATE_DEGREES)
        self.assertEqual(worst % 2, 0)
        self.assertEqual(
            worst // 2 * 45, int(raw["deficiency_budget_max"]["45"])
        )

    def test_build_states_covers_the_whole_deficiency_rectangle(self):
        edge_bounds = {order: (edges - 1, edges) for order, edges in (
            (20, 100), (21, 107), (22, 114), (23, 122), (24, 132))}
        catalog_windows = {
            (order, edges): (3 * order, 4 * order)
            for order, (low, high) in edge_bounds.items()
            for edges in (low, high)
        }
        states = build_states(unit_windows(), edge_bounds, catalog_windows)
        self.assertEqual(
            [(s.d, s.a, s.b) for s in states],
            [(d, a, b) for d in STATE_DEGREES
             for a in (0, 1) for b in (0, 1)],
        )
        self.assertEqual(states, sorted(states))
        for state in states:
            self.assertIsInstance(state, State)
            self.assertEqual(state.deficiency, state.a + state.b)
            self.assertEqual(
                state.excess_balance,
                2 * (state.a + state.b) - budget2(state.d, edge_bounds),
            )
            self.assertLessEqual(state.g_lo, state.g_hi)
            for value in dataclasses.astuple(state):
                self.assertIs(type(value), int)
        for d in STATE_DEGREES:
            zero = next(s for s in states if (s.d, s.a, s.b) == (d, 0, 0))
            self.assertEqual(zero.excess_balance, -budget2(d, edge_bounds))

    def test_build_states_keeps_pairs_whose_catalog_class_is_missing(self):
        windows = r35_edge_windows(DATA / "VALIDATION.json")
        edge_bounds = {order: (edges, edges) for order, edges in (
            (20, 100), (21, 107), (22, 114), (23, 122), (24, 132))}
        exact24 = [
            q_exact(adj)
            for adj in catalog_sample(DATA / "r45_24.g6", 132, 20)
        ]
        tight = {(24, 132): (min(exact24), max(exact24))}
        partial = {
            (order, bounds[1]): (3 * order, 4 * order)
            for order, bounds in edge_bounds.items() if order != 24
        }
        outer_states = build_states(windows, edge_bounds, partial)
        tight_states = build_states(windows, edge_bounds, partial | tight)
        self.assertEqual([s.d for s in outer_states], list(STATE_DEGREES))
        self.assertEqual([s.d for s in tight_states], list(STATE_DEGREES))
        for outer, inner in zip(outer_states, tight_states):
            self.assertEqual((outer.d, outer.a, outer.b),
                             (inner.d, inner.a, inner.b))
            self.assertLessEqual(outer.g_lo, inner.g_lo)
            self.assertLessEqual(inner.g_hi, outer.g_hi)
            if outer.d in (20, 24):  # the (24, 132) class is the missing one
                self.assertGreater(outer.g_hi - outer.g_lo,
                                   inner.g_hi - inner.g_lo)
            else:
                self.assertEqual((outer.g_lo, outer.g_hi),
                                 (inner.g_lo, inner.g_hi))

    def test_build_states_omits_only_provably_infeasible_edge_classes(self):
        edge_bounds = {20: (100, 100), 21: (107, 107), 22: (114, 114),
                       23: (122, 122), 24: (200, 200)}
        catalog_windows = {
            (order, bounds[1]): (3 * order, 4 * order)
            for order, bounds in edge_bounds.items() if order != 24
        }
        states = build_states(unit_windows(), edge_bounds, catalog_windows)
        # 2*200 > 24*13, so no integer histogram exists for (24, 200); the
        # d=20 and d=24 pairs are refuted, every other degree survives.
        self.assertEqual([s.d for s in states], [21, 22, 23])

    def test_state_interval_is_integral_on_the_outer_relaxation(self):
        windows = r35_edge_windows(DATA / "VALIDATION.json")
        edge_bounds = load_edge_bounds(DATA / "structural_tables.json")
        low, high = q_outer_interval(24, 132, windows)
        self.assertIs(type(low), int)
        self.assertIs(type(high), int)
        self.assertLessEqual(low, high)
        state = m3_state_interval(20, 0, 0, windows, edge_bounds)
        for value in dataclasses.astuple(state):
            self.assertIs(type(value), int)

    def test_infeasible_edge_class_is_rejected(self):
        windows = unit_windows()
        with self.assertRaises(InfeasibleEdgeClass) as caught:
            q_outer_interval(24, 200, windows)
        self.assertIn("order=24", str(caught.exception))
        self.assertIsInstance(caught.exception, ValueError)

    def test_missing_r35_window_is_rejected(self):
        with self.assertRaises(ValueError):
            q_outer_interval(20, 100, {k: (0, 0) for k in range(10)})

    def test_state_interval_rejects_out_of_range_arguments(self):
        windows = r35_edge_windows(DATA / "VALIDATION.json")
        edge_bounds = load_edge_bounds(DATA / "structural_tables.json")
        for d, a, b in ((20, -1, 0), (20, 0, -1), (20, 33, 0), (20, 0, 17),
                        (19, 0, 0), (25, 0, 0)):
            with self.assertRaises(ValueError):
                m3_state_interval(d, a, b, windows, edge_bounds)

    def test_load_edge_bounds_rejects_malformed_tables(self):
        table = json.loads(
            (DATA / "structural_tables.json").read_text()
        )["extremal_edge_tables_R45"]
        broken_min_max = {k: dict(v) for k, v in table.items()}
        broken_min_max["22"] = {"min": 115, "max": 114}
        fractional = {k: dict(v) for k, v in table.items()}
        fractional["22"] = {"min": 88, "max": 114.5}
        incomplete = {k: dict(v) for k, v in table.items() if k != "22"}
        cases = (
            {},
            {"extremal_edge_tables_R45": broken_min_max},
            {"extremal_edge_tables_R45": fractional},
            {"extremal_edge_tables_R45": incomplete},
            {"extremal_edge_tables_R45": {"22": {"min": 88}}},
        )
        with tempfile.TemporaryDirectory() as tmp:
            for index, payload in enumerate(cases):
                path = pathlib.Path(tmp) / f"bad{index}.json"
                path.write_text(json.dumps(payload), encoding="ascii")
                with self.assertRaises(ValueError):
                    load_edge_bounds(path)

    def test_r35_edge_windows_rejects_incomplete_validation(self):
        raw = json.loads((DATA / "VALIDATION.json").read_text())
        r35 = [r for r in raw["results"] if (r["s"], r["t"]) == (3, 5)]
        truncated = [r for r in r35 if int(next(iter(r["n_values"]))) < 12]
        unvalidated = [dict(r) for r in r35]
        unvalidated[-1]["ok"] = False
        conflicting = [dict(r) for r in r35] + [
            dict(r35[-1], edge_min=25, edge_max=25)
        ]
        cases = ({}, {"results": truncated}, {"results": unvalidated},
                 {"results": conflicting})
        with tempfile.TemporaryDirectory() as tmp:
            for index, payload in enumerate(cases):
                path = pathlib.Path(tmp) / f"bad{index}.json"
                path.write_text(json.dumps(payload), encoding="ascii")
                with self.assertRaises(ValueError):
                    r35_edge_windows(path)

    def test_q_outer_interval_exact_pin_both_intervals(self):
        # A five-vertex, five-edge class under deliberately trivial windows.
        windows = unit_windows()
        lo, hi = q_outer_interval(5, 5, windows)
        self.assertEqual((lo, hi), (4, 13))

        # Real R(3,5) windows; pin the order-24 boundary class.
        windows = r35_edge_windows(DATA / "VALIDATION.json")
        lo, hi = q_outer_interval(24, 132, windows)
        self.assertEqual((lo, hi), (1144, 1284))
        # The order-24 upper extremizer has lower triangle-sum 457,
        # so ceil(457/3)=153 and a floor mutation would widen q_hi to 1285.

    def test_build_states_raises_on_invalid_windows_and_catalog(self):
        edge_bounds = {order: (edges, edges) for order, edges in (
            (20, 100), (21, 107), (22, 114), (23, 122), (24, 132))}
        # Case (a): R(3,5) windows missing in-range degrees
        incomplete_windows = {k: (0, 0) for k in range(10)}
        catalog_windows = {
            (order, edges[1]): (3 * order, 4 * order)
            for order, edges in edge_bounds.items() if order != 20
        }
        with self.subTest("missing R(3,5) degree windows"):
            with self.assertRaises(ValueError):
                build_states(incomplete_windows, edge_bounds, catalog_windows)

        # Case (b): exact catalog window with non-int endpoint
        windows = unit_windows()
        bad_catalog_non_int = {(20, 100): (300, 300.5)}
        with self.subTest("non-integer catalog endpoint"):
            with self.assertRaises(ValueError):
                build_states(windows, edge_bounds, bad_catalog_non_int)

        # Case (c): exact catalog window with inverted endpoints
        bad_catalog_inverted = {(20, 100): (305, 304)}
        with self.subTest("inverted catalog window"):
            with self.assertRaises(ValueError):
                build_states(windows, edge_bounds, bad_catalog_inverted)


class CatalogQHistogramsTests(unittest.TestCase):
    def test_catalog_q_histograms_exact_two_classes_with_expected_edges_none(self):
        """5-cycle (q=5) and empty 5-vertex graph (q=0) split by observed edges."""
        cycle_adj = [0b10010, 0b00101, 0b01010, 0b10100, 0b01001]
        empty_adj = [0b00000, 0b00000, 0b00000, 0b00000, 0b00000]

        with tempfile.NamedTemporaryFile(mode='w', suffix='.g6', delete=False) as f:
            temp_path = f.name
            f.write(encode_graph6(5, cycle_adj) + '\n')
            f.write(encode_graph6(5, empty_adj) + '\n')

        try:
            result = catalog_q_histograms([(temp_path, 5, None, 2)])
            self.assertEqual(result[(5, 5)], {5: 1})
            self.assertEqual(result[(5, 0)], {0: 1})
        finally:
            pathlib.Path(temp_path).unlink(missing_ok=True)

    def test_catalog_q_histograms_invalid_inputs(self):
        """Invalid record shapes, empty catalog, count/order/edge mismatches, duplicates."""
        cycle_adj = [0b10010, 0b00101, 0b01010, 0b10100, 0b01001]

        with self.subTest("empty catalog"):
            with self.assertRaises(ValueError):
                catalog_q_histograms([])

        with self.subTest("malformed record (3-tuple)"):
            with self.assertRaises(ValueError):
                catalog_q_histograms([(None, 5, 5)])

        with tempfile.NamedTemporaryFile(mode='w', suffix='.g6', delete=False) as f:
            temp_path = f.name
            f.write(encode_graph6(5, cycle_adj) + '\n')

        try:
            with self.subTest("wrong expected count"):
                with self.assertRaises(ValueError):
                    catalog_q_histograms([(temp_path, 5, 5, 5)])

            with self.subTest("wrong expected order"):
                with self.assertRaises(ValueError):
                    catalog_q_histograms([(temp_path, 6, 5, 1)])

            with self.subTest("fixed-edge mismatch"):
                with self.assertRaises(ValueError):
                    catalog_q_histograms([(temp_path, 5, 3, 1)])
        finally:
            pathlib.Path(temp_path).unlink(missing_ok=True)

        with tempfile.NamedTemporaryFile(mode='w', suffix='.g6', delete=False) as f:
            temp_path1 = f.name
            f.write(encode_graph6(5, cycle_adj) + '\n')

        with tempfile.NamedTemporaryFile(mode='w', suffix='.g6', delete=False) as f:
            temp_path2 = f.name
            f.write(encode_graph6(5, cycle_adj) + '\n')

        try:
            with self.subTest("duplicate edge class"):
                with self.assertRaises(ValueError):
                    catalog_q_histograms([(temp_path1, 5, 5, 1), (temp_path2, 5, 5, 1)])
        finally:
            pathlib.Path(temp_path1).unlink(missing_ok=True)
            pathlib.Path(temp_path2).unlink(missing_ok=True)


if __name__ == "__main__":
    unittest.main()
