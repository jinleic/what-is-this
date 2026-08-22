# Higher-Order Identity $m=3$ Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use `subagent-driven-development` (recommended) or `executing-plans` to implement this plan task-by-task. Steps use checkbox (`- [x]`) syntax for tracking.

**Goal:** Determine, with exact arithmetic and an independent checker, whether the McKay–Radziszowski $m=3$ neighborhood identity proves a publishable positive-deficiency restriction on every hypothetical Ramsey$(5,5,45)$ graph.

**Architecture:** First implement and exhaustively validate the universal $m=3$ identity. Then compute exact $q=\sum_x\binom{d_x}{2}-t$ histograms for every published-complete $R(4,5,r,e)$ catalog at orders 17–24; use those exact windows where available and a sound degree-histogram outer envelope from the complete $R(3,5,k)$ catalogs everywhere else. Build the finite $n=45$ local-state cone and search for exact affine certificates of the approved consequences. A separate checker reconstructs the cone without importing the producer. Solver output proposes coefficients only; exact rational inequalities decide acceptance.

**Tech Stack:** Python 3.14.3, standard-library `unittest`, integer bitsets, `fractions.Fraction`, SciPy 1.18.0 `linprog` for coefficient discovery only, validated graph6 catalogs.

## Global Constraints

- Canonical design: `math/r55/bench/spec.json:next_research_campaign`.
- AI and SciPy remain outside the trust base; every accepted coefficient is rationalized and checked exactly.
- Reuse `math/r55/src/check_ramsey.py` for graph6 parsing and `popcount`; do not create a second parser.
- Missing inputs, count drift, identity mismatch, or certificate mismatch must exit nonzero.
- Do not change the frozen benchmark engine or candidate-3 artifacts.
- Do not add dependencies, run raw catalog generation, or start the deferred SDP/pair-slack tracks.
- No commit, push, or external mutation without separate authorization; this project tree is currently untracked by Git.
- A valid negative result is `M3_NO_CUT_IN_FROZEN_CONE`, backed by exact primal witnesses for every instantiable route; never relabel a numerical trend, unresolved certification, or relaxation result as a theorem about all possible uses of $m=3$.

---

### Task 1: Exact $m=3$ identity kernel

**Files:**
- Create: `math/r55/src/subgraph_identities.py`
- Create: `math/r55/tests/test_subgraph_identities.py`

**Interfaces:**
- Consumes: `check_ramsey.popcount`.
- Produces:
  - `induced_subgraph(adj: list[int], vertices: list[int]) -> list[int]`
  - `complement(adj: list[int]) -> list[int]`
  - `edge_count(adj: list[int]) -> int`
  - `triangle_count(adj: list[int]) -> int`
  - `induced_p3_count(adj: list[int]) -> int`
  - `independent_triple_count(adj: list[int]) -> int`
  - `motif3(adj: list[int]) -> Motif3`
  - `m3_vertex_scaled(n: int, neighborhood: list[int], dual: list[int]) -> int`
  - `m3_residual_scaled(adj: list[int]) -> int`

The exact identity is stored after multiplication by three:

$$
3\sum_v t(G_v^-)=\sum_v\left[(n+3-3d_v)e(G_v^+)+3p_3(G_v^+)+6t(G_v^+)\right].
$$

The module docstring must include the independent counting derivation after
dividing the displayed identity by three. In that unscaled form, the left side
counts induced $K_3\sqcup K_1$ subgraphs. Grouping the first term on the right
by a triangle $abc$ gives
$\sum_{\triangle abc}\sum_{x\notin\{a,b,c\}}(1-r_x)$, where $r_x$ is the
number of neighbors of $x$ in that triangle. Its induced four-vertex
contributions are $+1$ for $K_3\sqcup K_1$, $0$ for a paw, $-2$ for a
diamond, and $-8$ for $K_4$. The sums of $p_3(G_v^+)$ and
$2t(G_v^+)$ contribute $+2$ per diamond and $+8$ per $K_4$, respectively,
leaving exactly the unscaled left side. Multiplying back by three yields the
stored integer identity. The exhaustive labeled-graph test is a second,
mechanically independent check of this derivation.

- [x] **Step 1: Write the exhaustive failing tests**

Create `test_subgraph_identities.py` with a local unlabeled-graph constructor and these observable contracts:

```python
import itertools
import pathlib
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from subgraph_identities import (  # expected to fail before Task 1 implementation
    complement, independent_triple_count, induced_p3_count, induced_subgraph,
    m3_residual_scaled, motif3, triangle_count,
)


def graph_from_mask(n, mask):
    adj = [0] * n
    for bit, (u, v) in enumerate(itertools.combinations(range(n), 2)):
        if (mask >> bit) & 1:
            adj[u] |= 1 << v
            adj[v] |= 1 << u
    return adj


class M3IdentityTests(unittest.TestCase):
    def test_motif_counts(self):
        # Build C5 explicitly so graph-bit ordering cannot hide a test error.
        c5 = [0] * 5
        for u, v in ((0, 1), (1, 2), (2, 3), (3, 4), (4, 0)):
            c5[u] |= 1 << v
            c5[v] |= 1 << u
        self.assertEqual(motif3(c5), (5, 5, 0, 5, 0))
        empty4 = [0] * 4
        self.assertEqual(independent_triple_count(empty4), 4)

    def test_complement_is_involution(self):
        for n in range(1, 6):
            for mask in range(1 << (n * (n - 1) // 2)):
                adj = graph_from_mask(n, mask)
                self.assertEqual(complement(complement(adj)), adj)

    def test_m3_identity_exhaustive_through_six_vertices(self):
        for n in range(1, 7):
            for mask in range(1 << (n * (n - 1) // 2)):
                self.assertEqual(m3_residual_scaled(graph_from_mask(n, mask)), 0)

    def test_mutated_triangle_coefficient_is_detected(self):
        k4 = graph_from_mask(4, (1 << 6) - 1)
        self.assertEqual(m3_residual_scaled(k4), 0)
        # The production identity uses coefficient 6. K4 detects replacing it by 5.
        n = 4
        mutated = 0
        for v in range(n):
            neighborhood_vertices = [
                u for u in range(n) if (k4[v] >> u) & 1
            ]
            dual_vertices = [
                u for u in range(n)
                if u != v and not ((k4[v] >> u) & 1)
            ]
            neighborhood = induced_subgraph(k4, neighborhood_vertices)
            dual = induced_subgraph(k4, dual_vertices)
            f = motif3(neighborhood)
            mutated += 3 * triangle_count(dual) - (
                (n + 3 - 3 * len(neighborhood)) * f.edges
                + 3 * f.induced_p3 + 5 * f.triangles
            )
        self.assertNotEqual(mutated, 0)


if __name__ == "__main__":
    unittest.main()
```


- [x] **Step 2: Run the test and observe the intended failure**

Run:

```bash
cd /Users/jinleic/jinleic-workspace/math
./.venv/bin/python -m unittest r55/tests/test_subgraph_identities.py -v
```

Expected: import failure for `subgraph_identities`.

- [x] **Step 3: Implement the exact kernel**

Use this implementation shape in `subgraph_identities.py`:

```python
#!/usr/bin/env python3
from collections import namedtuple

from check_ramsey import popcount

Motif3 = namedtuple(
    "Motif3", "order edges triangles induced_p3 independent_triples"
)


def complement(adj):
    n = len(adj)
    full = (1 << n) - 1
    return [full & ~adj[v] & ~(1 << v) for v in range(n)]


def induced_subgraph(adj, vertices):
    out = [0] * len(vertices)
    for i, u in enumerate(vertices):
        for j in range(i + 1, len(vertices)):
            v = vertices[j]
            if (adj[u] >> v) & 1:
                out[i] |= 1 << j
                out[j] |= 1 << i
    return out


def edge_count(adj):
    return sum(popcount(row) for row in adj) // 2


def triangle_count(adj):
    total = 0
    for u in range(len(adj)):
        later = adj[u] & ~((1 << (u + 1)) - 1)
        while later:
            bit = later & -later
            v = bit.bit_length() - 1
            total += popcount(adj[u] & adj[v])
            later ^= bit
    return total // 3


def induced_p3_count(adj):
    triangles = triangle_count(adj)
    wedges = sum(popcount(row) * (popcount(row) - 1) // 2 for row in adj)
    return wedges - 3 * triangles


def independent_triple_count(adj):
    return triangle_count(complement(adj))


def motif3(adj):
    return Motif3(
        len(adj), edge_count(adj), triangle_count(adj), induced_p3_count(adj),
        independent_triple_count(adj),
    )


def m3_vertex_scaled(n, neighborhood, dual):
    if len(neighborhood) + len(dual) != n - 1:
        raise ValueError("neighborhood and dual orders must sum to n-1")
    d = len(neighborhood)
    f = motif3(neighborhood)
    return 3 * triangle_count(dual) - (
        (n + 3 - 3 * d) * f.edges + 3 * f.induced_p3 + 6 * f.triangles
    )


def m3_residual_scaled(adj):
    n = len(adj)
    return sum(
        m3_vertex_scaled(
            n,
            induced_subgraph(adj, [u for u in range(n) if (adj[v] >> u) & 1]),
            induced_subgraph(
                adj, [u for u in range(n) if u != v and not ((adj[v] >> u) & 1)]
            ),
        )
        for v in range(n)
    )
```

- [x] **Step 4: Run exhaustive tests and the existing graph-parser selftest**

Run:

```bash
cd /Users/jinleic/jinleic-workspace/math
./.venv/bin/python -m unittest r55/tests/test_subgraph_identities.py -v
./.venv/bin/python r55/src/check_ramsey.py --selftest
```

Expected: all new tests `OK`; existing selftest reports `selftest OK`.

---

### Task 2: Sound motif envelopes and finite $n=45$ state cone

**Files:**
- Create: `math/r55/src/m3_deficiency_cone.py`
- Create: `math/r55/tests/test_m3_deficiency_cone.py`
- Create at execution time: `math/r55/data/higher_identity_m3.json`

**Interfaces:**
- Consumes: Task 1 primitives, `data/structural_tables.json`, `data/VALIDATION.json`, `data/VALIDATION_extreme.json`, complete `r35_*.g6` catalogs, every validated $R(4,5)$ edge-extremal catalog at orders 17–23, and the complete `r45_24.g6` catalog.
- Produces:
  - `load_edge_bounds(path) -> dict[int, tuple[int, int]]`
  - `r35_edge_windows(validation_path) -> dict[int, tuple[int, int]]`
  - `degree_histograms(order: int, edges: int)` returns integer-count tuples
  - `q_outer_interval(order: int, edges: int, windows) -> tuple[int, int]`, where $q=\sum_x\binom{d_x}{2}-t$
  - `catalog_q_histograms(catalog_records) -> dict[tuple[int, int], dict[int, int]]`
  - `q_interval(order, edges, windows, catalog_windows) -> tuple[int, int]`
  - `m3_state_interval(d, a, b, windows, edge_bounds, catalog_windows=None) -> State`
  - `build_states(windows, edge_bounds, catalog_windows) -> list[State]`
  - deterministic analysis JSON with source hashes, 656-graph replay, exact per-catalog $q$ histograms, the published-completeness dependency, and every state interval.

For an $R(4,5,r,e)$ graph $H$:

$$
i_3(H)=\binom r3-e(r-2)+\sum_x\binom{d_x}{2}-t(H).
$$

Each $N_H(x)$ is an $R(3,5,d_x)$ graph, so complete $R(3,5,k)$ edge windows give

$$
\sum_x e_{35}^{\min}(d_x)\le 3t(H)\le
\sum_x e_{35}^{\max}(d_x).
$$

Enumerating all graphical *degree histograms is not required for soundness*; enumerate every integer histogram satisfying the Ramsey degree window and handshake equation. This is an outer relaxation, hence safe. Do not silently impose graphicality.

For a catalog explicitly claimed complete by the Angeltveit–McKay data source,
the exact minimum and maximum observed $q$ replace the outer interval. For
missing edge classes, the outer interval remains authoritative. This tightens
the cone without assuming that an unlisted positive-deficiency class is empty.
The published completeness of the order-24 and edge-extremal catalogs remains
a named trust-root dependency; local hashing and validation do not re-prove it.

- [x] **Step 1: Write failing envelope tests**

Tests must cover:

```python
import itertools
import pathlib
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from check_ramsey import parse_graph6_line
from m3_deficiency_cone import (
    degree_histograms, load_edge_bounds, m3_state_interval, q_interval,
    q_outer_interval, r35_edge_windows,
)
from subgraph_identities import complement, edge_count, m3_residual_scaled, motif3


def read_graphs(path, limit=None, required_edges=None):
    graphs = []
    with path.open(encoding="ascii") as handle:
        for line in handle:
            if not line.strip():
                continue
            _, adj = parse_graph6_line(line)
            if required_edges is not None and edge_count(adj) != required_edges:
                continue
            graphs.append(adj)
            if limit is not None and len(graphs) == limit:
                break
    return graphs


def q_exact(adj):
    f = motif3(adj)
    return sum(
        row.bit_count() * (row.bit_count() - 1) // 2 for row in adj
    ) - f.triangles


class M3ConeTests(unittest.TestCase):
    def test_degree_histograms_obey_order_and_handshake(self):
        for hist in degree_histograms(5, 5):
            lo = max(0, 5 - 18)
            self.assertEqual(sum(hist), 5)
            self.assertEqual(sum((lo + i) * c for i, c in enumerate(hist)), 10)

    def test_outer_interval_contains_exact_c5_value(self):
        windows = {k: (0, k * (k - 1) // 2) for k in range(14)}
        lo, hi = q_outer_interval(5, 5, windows)
        exact_q = 5  # five degree-2 wedges, zero triangles
        self.assertLessEqual(lo, exact_q)
        self.assertLessEqual(exact_q, hi)

    def test_exact_catalog_window_overrides_outer_relaxation(self):
        windows = {k: (0, k * (k - 1) // 2) for k in range(14)}
        catalog_windows = {(5, 5): (7, 7)}
        self.assertEqual(q_interval(5, 5, windows, catalog_windows), (7, 7))

    def test_state_interval_contains_exact_catalog_pairs(self):
        windows = r35_edge_windows(ROOT / "data" / "VALIDATION.json")
        edge_bounds = load_edge_bounds(ROOT / "data" / "structural_tables.json")
        max_edges = {order: bounds[1] for order, bounds in edge_bounds.items()}
        pairs = (
            (20, 100, ROOT / "data" / "r45extreme" / "r4520.100.g6",
             24, 132, ROOT / "data" / "r45_24.g6"),
            (21, 107, ROOT / "data" / "r45extreme" / "r4521.107.g6",
             23, 122, ROOT / "data" / "r45extreme" / "r4523.122.g6"),
            (22, 114, ROOT / "data" / "r45extreme" / "r4522.114.g6",
             22, 114, ROOT / "data" / "r45extreme" / "r4522.114.g6"),
        )
        for d, ex, x_path, m, ey, y_path in pairs:
            self.assertEqual(m, 44 - d)
            state = m3_state_interval(
                d, max_edges[d] - ex, max_edges[m] - ey, windows, edge_bounds
            )
            for x, y in itertools.product(
                read_graphs(x_path, 20, ex), read_graphs(y_path, 20, ey)
            ):
                exact = (
                    3 * (m * (m - 1) * (m - 2) // 6 - ey * (m - 2) + q_exact(y))
                    - (48 - 3 * d) * ex - 3 * q_exact(x)
                )
                self.assertLessEqual(state.g_lo, exact)
                self.assertLessEqual(exact, state.g_hi)

    def test_all_656_known_r55_graphs_have_zero_m3_residual(self):
        published = read_graphs(ROOT / "data" / "r55_42some.g6")
        self.assertEqual(len(published), 328)
        tested = 0
        for adj in published:
            self.assertEqual(m3_residual_scaled(adj), 0)
            self.assertEqual(m3_residual_scaled(complement(adj)), 0)
            tested += 2
        self.assertEqual(tested, 656)
```

- [x] **Step 2: Run tests and observe missing-interface failures**

Run:

```bash
cd /Users/jinleic/jinleic-workspace/math
./.venv/bin/python -m unittest r55/tests/test_m3_deficiency_cone.py -v
```

Expected: import failures for the cone interfaces.

- [x] **Step 3: Implement degree-histogram and motif-envelope logic**

Use exact integers. The core must implement these formulas:

```python
import json
from dataclasses import dataclass


@dataclass(frozen=True, order=True)
class State:
    d: int
    a: int
    b: int
    deficiency: int
    excess_balance: int
    g_lo: int
    g_hi: int


def load_edge_bounds(path):
    raw = json.loads(path.read_text())
    return {
        int(order): (int(values["min"]), int(values["max"]))
        for order, values in raw["extremal_edge_tables_R45"].items()
    }


def budget2(d, edge_bounds):
    m = 44 - d
    max_edges = {order: bounds[1] for order, bounds in edge_bounds.items()}
    h2_min = (
        2 * (m * (m - 1) // 2 - max_edges[m] - max_edges[d])
        - d * (45 - 2 * d)
    )
    return -h2_min


def degree_bounds(order):
    return max(0, order - 18), min(13, order - 1)


def q_outer_interval(order, edges, windows):
    q_lo = None
    q_hi = None
    for hist in degree_histograms(order, edges):
        first = degree_bounds(order)[0]
        wedge_sum = sum(
            count * degree * (degree - 1) // 2
            for degree, count in enumerate(hist, start=first)
        )
        tri3_lo = sum(count * windows[degree][0]
                      for degree, count in enumerate(hist, start=first))
        tri3_hi = sum(count * windows[degree][1]
                      for degree, count in enumerate(hist, start=first))
        triangle_lo = (tri3_lo + 2) // 3
        triangle_hi = tri3_hi // 3
        if triangle_lo > triangle_hi:
            continue
        local_lo = wedge_sum - triangle_hi
        local_hi = wedge_sum - triangle_lo
        q_lo = local_lo if q_lo is None else min(q_lo, local_lo)
        q_hi = local_hi if q_hi is None else max(q_hi, local_hi)
    if q_lo is None:
        raise ValueError(f"no degree histogram for order={order}, edges={edges}")
    return q_lo, q_hi


def q_interval(order, edges, windows, catalog_windows):
    key = (order, edges)
    if key in catalog_windows:
        return catalog_windows[key]
    return q_outer_interval(order, edges, windows)


def m3_state_interval(
    d, a, b, windows, edge_bounds, catalog_windows=None
):
    catalog_windows = {} if catalog_windows is None else catalog_windows
    max_edges = {order: bounds[1] for order, bounds in edge_bounds.items()}
    m = 44 - d
    ex = max_edges[d] - a
    ey = max_edges[m] - b
    qx_lo, qx_hi = q_interval(d, ex, windows, catalog_windows)
    qy_lo, qy_hi = q_interval(m, ey, windows, catalog_windows)
    dual_base = 3 * (m * (m - 1) * (m - 2) // 6 - ey * (m - 2))
    rhs_base = (48 - 3 * d) * ex
    g_lo = dual_base + 3 * qy_lo - rhs_base - 3 * qx_hi
    g_hi = dual_base + 3 * qy_hi - rhs_base - 3 * qx_lo
    return State(
        d=d, a=a, b=b, deficiency=a + b,
        excess_balance=2 * (a + b) - budget2(d, edge_bounds),
        g_lo=g_lo, g_hi=g_hi,
    )
```

Here $K=\overline{G_v^-}$, $e(G_v^+)=E_{45}(d)-a$, and
$e(K)=E_{45}(m)-b$. Substitution into the existing $m=2$ excess identity
gives `excess_balance = 2(a+b)-budget2(d)` and therefore
$\sum_v excess\_balance(v)=0$. Task 1 gives an actual $m=3$ contribution
$g_v\in[g_{lo},g_{hi}]$ with $\sum_v g_v=0$, hence
$\sum_v g_{lo}\le0\le\sum_v g_{hi}$. These are the only global equalities and
inequalities used by the affine certificate.

`degree_histograms` must enumerate counts for every degree in the inclusive Ramsey window, prune only by remaining vertex count and remaining degree sum, and yield in lexicographic order.

`r35_edge_windows` must synthesize order zero as `(0, 0)`. `build_states` must
enumerate every $d\in\{20,\ldots,24\}$,
$0\le a\le E_{45}(d)-e_{45}(d)$, and
$0\le b\le E_{45}(44-d)-e_{45}(44-d)$. It may omit a pair only when the
outer degree-histogram constraints prove one side infeasible. It must never
treat a missing $R(4,5)$ catalog as infeasibility. Edge extrema have one
semantic owner, `data/structural_tables.json`; the new module loads them and
derives `budget2` rather than copying constants.

- [x] **Step 4: Add fail-closed input provenance and catalog calibration**

The CLI must:

1. hash every input file with SHA-256;
2. load edge windows from the validated complete `r35_k.g6` records and synthesize order zero;
3. assert all 328 lines plus complements give 656 zero residuals;
4. stream every validated edge-extremal file at orders 17–23 and every graph in the complete order-24 catalog;
5. persist the exact $q$ frequency histogram for each `(order, edges)` class;
6. verify every exact $q(H)$ lies inside `q_outer_interval(order,e)`;
7. use exact catalog windows for present complete classes and outer windows for missing classes;
8. write JSON only after all checks pass, using sorted keys and deterministic list ordering.

The JSON schema is fixed before search:

```text
schema_version: 1
campaign_id: higher_order_identity_positive_deficiency_m3
inputs: sorted records {relative_path, sha256, bytes, graph_count}
trust_roots: sorted published-completeness statements and source URLs
replay: {published_graphs: 328, complements: 328, residual_zero: 656}
r35_edge_windows: sorted records {order, edge_min, edge_max}
catalog_q_histograms: sorted records
  {order, edges, source, graph_count, q_counts: sorted [{q, count}]}
catalog_windows: sorted records {order, edges, q_min, q_max, source}
states: sorted records
  {d, a, b, deficiency, excess_balance, g_lo, g_hi,
   x_interval_source, y_interval_source}
searches: [] until Task 3
disposition: M3_ANALYZED until Task 3
```

No wall-clock timestamp or absolute path enters the canonical JSON. The CLI
writes to a temporary sibling, parses it back, then atomically replaces the
target only after schema and count checks pass.

- [x] **Step 5: Run focused tests and full analysis**

Run:

```bash
cd /Users/jinleic/jinleic-workspace/math
./.venv/bin/python -m unittest r55/tests/test_m3_deficiency_cone.py -v
./.venv/bin/python r55/src/m3_deficiency_cone.py \
  --output r55/data/higher_identity_m3.json
```

Expected terminal contract:

```text
M3 REPLAY: 656/656 residual zero
M3 ENVELOPES: catalog violations=0
M3 STATES: N states (N >= 1)
WROTE r55/data/higher_identity_m3.json
```

---

### Task 3: Exact affine cuts and obstruction witnesses

**Files:**
- Create: `math/r55/src/search_m3_cuts.py`
- Create: `math/r55/tests/test_search_m3_cuts.py`
- Modify at execution time: `math/r55/data/higher_identity_m3.json`

**Interfaces:**
- Consumes: Task 2 states and a frozen objective registry.
- Produces an immutable `Certificate` with fields `objective_id`, `alpha`, `beta`, `gamma`, and `delta`; rational values serialize as reduced fraction strings.
- Produces an immutable `PrimalWitness` with `objective_id` and sorted `(state_index, weight)` rational pairs.
- Produces `exact_affine_bound(states, objective_id) -> Certificate`, `verify_certificate(states, certificate)`, `verify_primal_witness(states, witness)`, and `verify_search_record(states, record)`.
- Every certificate satisfies $\gamma,\delta\ge0$ and, for every state,

$$
objective(state)\le \alpha+\beta\,balance(state)
 +\gamma g_{lo}(state)-\delta g_{hi}(state).
$$

Summing over 45 vertices uses $\sum balance=0$, $\sum g_{lo}\le0$, and $\sum g_{hi}\ge0$, proving the exact global bound $\sum objective\le45\alpha$.

- [x] **Step 1: Write failing synthetic-certificate tests**

Cover:

```python
import pathlib
import sys
import unittest
from dataclasses import replace
from fractions import Fraction

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from m3_deficiency_cone import State
from search_m3_cuts import (
    PrimalWitness, exact_affine_bound, verify_certificate,
    verify_primal_witness,
)


class ExactCutTests(unittest.TestCase):
    def test_rationalized_certificate_checks_every_state(self):
        states = [
            State(20, 0, 0, 0, -12, -4, 2),
            State(24, 0, 0, 0, -12, -2, 4),
        ]
        cert = exact_affine_bound(states, "total_deficiency")
        verify_certificate(states, cert)

    def test_tampered_alpha_is_rejected(self):
        states = [
            State(20, 0, 0, 0, -12, -4, 2),
            State(24, 0, 0, 0, -12, -2, 4),
        ]
        cert = exact_affine_bound(states, "total_deficiency")
        tampered = replace(cert, alpha=str(Fraction(cert.alpha) - 1))
        with self.assertRaises(ValueError):
            verify_certificate(states, tampered)

    def test_exact_primal_rejection_witness(self):
        states = [State(20, 0, 0, 0, 0, -1, 1)]
        witness = PrimalWitness("degree20_count", ((0, "45"),))
        self.assertEqual(verify_primal_witness(states, witness), Fraction(45))

    def test_wrong_objective_id_is_rejected(self):
        states = [State(20, 0, 0, 0, 0, -1, 1)]
        witness = PrimalWitness("not_registered", ((0, "45"),))
        with self.assertRaises(ValueError):
            verify_primal_witness(states, witness)
```

- [x] **Step 2: Run tests and observe missing-interface failures**

Run:

```bash
cd /Users/jinleic/jinleic-workspace/math
./.venv/bin/python -m unittest r55/tests/test_search_m3_cuts.py -v
```

- [x] **Step 3: Implement exact upper certificates and lower witnesses**

Freeze this registry in source before the first full run:

```text
total_deficiency: objective=state.deficiency, accept upper_bound <= 315
degree20_count: objective=int(state.d == 20), accept upper_bound < 1
deficiency_ge8_count: objective=int(state.deficiency >= 8), accept upper_bound < 1
```

For an actual 45-vertex graph the two count objectives are nonnegative
integers, so an exact upper bound below one proves the count is zero. No
integrality assumption is used for the continuous relaxation itself.

Solve the continuous dual discovery LP with variables
$(\alpha,\beta,\gamma,\delta)$, bounds
`[(None, None), (None, None), (0, None), (0, None)]`, objective
`[1, 0, 0, 0]`, and one row
`[-1, -balance, -g_lo, g_hi] <= -objective(state)` per state. Convert the
returned $(\beta,\gamma,\delta)$ to `Fraction` using the frozen denominator
sequence `(100, 1000, 10000, 100000, 1000000)`. For each candidate, derive
$\alpha$ exactly:

```python
def exact_alpha(states, objective, beta, gamma, delta):
    return max(
        Fraction(objective(state))
        - beta * state.excess_balance
        - gamma * state.g_lo
        + delta * state.g_hi
        for state in states
    )
```

`verify_certificate` must resolve `objective_id` through the frozen registry,
reject unknown IDs, check coefficient signs and every state inequality, and
return the computed exact bound $45\alpha$. No caller-provided lambda or stored
boolean controls acceptance.

For every route not accepted by an exact upper certificate, solve the primal
discovery LP with variables $x_s\ge0$:

```text
maximize sum objective(s) * x_s
subject to sum x_s = 45
           sum excess_balance(s) * x_s = 0
           sum g_lo(s) * x_s <= 0
           sum g_hi(s) * x_s >= 0
```

Use HiGHS dual simplex only to identify a sparse basic support. First try the
indices with primal weight greater than `1e-9`. If exact reconstruction fails,
form one deterministic candidate set from those indices and the 32 states with
smallest `(absolute_dual_slack, state_index)`. Enumerate its subsets of at most
four states and every active/inactive choice for the two $g$ inequalities;
solve each basis by `Fraction` Gaussian elimination. Do not expand this set
after seeing results: failure becomes `CERTIFICATION_UNRESOLVED`.
`verify_primal_witness` independently requires nonnegative weights, total
weight 45, zero excess balance, aggregate $g_{lo}\le0\le g_{hi}$, and
recomputes the exact objective.

A route is `REJECTED_BY_EXACT_WITNESS` only when the witness value is $>315$
for `total_deficiency` or $\ge1$ for either count objective. Otherwise it is
`CERTIFICATION_UNRESOLVED`; numerical solver status or a rounded coefficient
never proves a negative result.

Search the three frozen objectives above. A fourth named `(d,a,b)` family is
searched only when a hash-pinned complete-gluing-cover manifest names it as
required. It is accepted only if an exact upper certificate gives
$45\alpha<1$ and the independent checker validates that necessity certificate.

The current repository has no complete $n=45$ gluing-cover manifest, so route 4
must be recorded as `UNAVAILABLE_NO_COVER_CERTIFICATE` unless such an
independently validated input is discovered before the frozen search begins.
Eliminating an arbitrary local family is not significant success. Do not tune
objective definitions, thresholds, denominator sequence, or support-expansion
rules after seeing results.

- [x] **Step 4: Run the cut search and persist exact decisions**

Run:

```bash
cd /Users/jinleic/jinleic-workspace/math
./.venv/bin/python r55/src/search_m3_cuts.py \
  --analysis r55/data/higher_identity_m3.json
```

Expected terminal contract is exactly one of:

```text
M3 CUT SEARCH COMPLETE status=ACCEPTED_CUT accepted=N (N >= 1)
```

```text
M3 CUT SEARCH COMPLETE status=NO_CUT_IN_FROZEN_CONE accepted=0
```

```text
M3 CUT SEARCH COMPLETE status=CERTIFICATION_UNRESOLVED accepted=0
```

`NO_CUT_IN_FROZEN_CONE` requires an exact rejection witness for every
instantiable route. `CERTIFICATION_UNRESOLVED` exits nonzero and cannot trigger
a scientific disposition. Every persisted search record contains its objective
ID, exact upper certificate and bound when available, exact lower witness and
value when available, derived route status, and no unchecked acceptance field.

---

### Task 4: Disjoint checker and evidence-grade disposition

**Files:**
- Create: `math/r55/src/check_m3_certificate.py`
- Create: `math/r55/tests/test_check_m3_certificate.py`
- Modify: `math/r55/bench/spec.json`
- Modify: `math/PROGRESS.md`
- Modify: `math/r55/README.md`
- Modify: `math/r55/notes/structural_constraints_2026-08-15.md`

**Interfaces:**
- Consumes: frozen validation files, graph catalogs, and `data/higher_identity_m3.json`.
- Produces an independent pass/fail verdict. It must not import `subgraph_identities`, `m3_deficiency_cone`, or `search_m3_cuts`.

- [x] **Step 1: Write checker corruption tests**

Create a temporary copy of the analysis artifact and independently test rejection of:

1. one changed source SHA-256;
2. one changed state endpoint;
3. one reduced `alpha` coefficient;
4. one changed primal weight;
5. one false exact bound or route status;
6. one changed 656-graph replay count.

The untouched artifact must pass. Separately exhaust all labeled graphs through
five vertices against the checker's own residual implementation and show that
changing the scaled triangle coefficient from 6 to 5 fails on $K_4$.

- [x] **Step 2: Implement the independent checker**

The checker deliberately reimplements the following small trust kernel:

- graph complement, induced subgraph, edge, triangle, and induced-$P_3$ counts;
- the scaled $m=3$ residual, recomputed as zero on all 328 published graphs and their complements;
- degree histogram enumeration;
- $R(3,5,k)$ edge-window extraction from `VALIDATION.json`;
- exact $q$ histograms from all applicable published-complete $R(4,5)$ files;
- $q=W-t$ outer interval arithmetic for every missing edge class;
- all `(d,a,b)` state endpoints and aggregate balance signs;
- the frozen objective registry, exact dual certificates, exact primal witnesses, bounds, thresholds, and derived route/disposition statuses;
- source hashes, catalog counts, and frozen replay counts.

It may import only Python standard-library modules and
`check_ramsey.parse_graph6_line`/`popcount`. It must not trust any producer
count, bound, boolean, or disposition without recomputation. Any mismatch
prints the first exact counterexample and exits 1.

- [x] **Step 3: Run all focused tests and the independent checker**

Run:

```bash
cd /Users/jinleic/jinleic-workspace/math
./.venv/bin/python -m unittest discover -s r55/tests -p 'test_*m3*.py' -v
./.venv/bin/python r55/src/check_m3_certificate.py \
  r55/data/higher_identity_m3.json
```

Expected: focused tests `OK`; checker prints `M3 EVIDENCE VERIFIED` and the exact disposition.

- [x] **Step 4: Apply the frozen scientific gate**

If at least one exact upper certificate meets a quantified success route:

- set `next_research_campaign.status` to `M3_ACCEPTED_PENDING_NOVELTY_AUDIT`;
- record the theorem statement with exact coefficients and hypotheses in the existing structural note;
- conduct the primary-source novelty audit before using “new” or “publishable”;
- do not claim $R(5,5)\le45$.

If every instantiable route is rejected by an exact primal witness:

- set status to `M3_NO_CUT_IN_FROZEN_CONE`;
- record each exact witness and state only that the frozen continuous cone cannot prove an accepted route;
- preserve all artifacts and hashes;
- transition to a separately reviewed $m=4$ plan, as required by the approved design.

If any route is `CERTIFICATION_UNRESOLVED`:

- set status to `M3_CERTIFICATION_UNRESOLVED`;
- preserve the verified upper/lower bounds without a positive or negative claim;
- stop before $m=4$ and record the exact unresolved gap and basis-reconstruction failure.

- [x] **Step 5: Final integrity verification**

Parse both JSON files, verify every recorded SHA-256, rerun the checker, and confirm no frozen benchmark hash or candidate artifact changed. Report the mathematical result first, then files and verification evidence.

---

## Plan Self-Review

- **Spec coverage:** Implements the approved first $m=3$ slice, independent derivation, exact arithmetic, 656-graph replay, exact catalog motif histograms plus sound missing-class envelopes, finite cone, exact upper certificates and lower witnesses, solver outside the trust base, an independent checker, all currently instantiable quantified success routes, and fail-closed positive, negative, or unresolved dispositions.
- **Placeholder scan:** No unresolved implementation placeholders, deferred code bodies, or production stubs remain in the plan.
- **Type consistency:** `State`, motif interfaces, state endpoints, certificate coefficient names, and artifact path are consistent across all four tasks.
- **Scope:** No $m=4$, mixed-identity, SDP, gluing-engine, benchmark, or raw-generation work enters this plan.
