# Higher-Order Identity $m=4$ Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use `subagent-driven-development` (recommended) or `executing-plans` to implement this plan task-by-task. Steps use checkbox (`- [x]`) syntax for tracking.

**Goal:** Determine, with exact arithmetic and an independent checker, whether the McKay–Radziszowski $m=4$ neighborhood identity proves a publishable positive-deficiency or degree restriction on every hypothetical Ramsey$(5,5,45)$ graph. A decisive negative (`M4_NO_CUT_IN_FROZEN_BASIS`) is an allowed, evidence-grade outcome that closes the $T$-family identity line at $n=45$.

**Motivation and strategic context.**

Literature status (DS1 rev Apr 24 2026): $R(5,5)\in\{43,\dots,46\}$ — lower bound 43 (Exoo 1989), upper bound $46$ (Angeltveit–McKay, arXiv:2409.15709 v2). At $n=45$ no published structural restriction exists beyond the degree window $\{20,\dots,24\}$; $n=46$ is governed by budget-plus-gluing lemmas; $n=47,48,49$ are eliminated, and it was the $m=4$ identity that killed $n=49$ in McKay–Radziszowski Thm 3.2 (1997), by hand.

The research agent verified the closed form, for graphs whose every neighborhood is $K_4$-free (in particular every $(5,5)$-graph, since a $K_4$ in $N(v)$ would complete a $K_5$ with $v$):

$$
12\sum_v i_4\big(\overline{D_v}\big)
=\sum_v\Big[3(n+2-2d_v)\,t(N_v)+4\,\mathrm{diamond}(N_v)-6\,s(K_3{+}K_1,N_v)\Big],
$$

by brute force on all labeled graphs $n\le 6$ plus random $7$–$9$. Plan drafting then re-derived and re-verified every constant independently in-repo (see Task 1 fixtures): the identity holds universally in the corrected form with a $+12\sum_v K_4(N_v)$ term on the right (§ Task 1), and the $n=49$ ground truth replays exactly on this repository's own frozen catalogs. Two published facts bound the algebraic line: $m\ge 5$ of the $T$-family is void for $(5,5)$-graphs, and MR97 Lemma 2.3 shows no further separable identities of degree $\le 6$. The next machinery tier is Engström's mixed-quadratic identity (arXiv:1002.4304), which this plan explicitly excludes.

Strategic scenario. The recognized bottleneck for $R(5,5)\le 45$ is the pointed-gluing enumeration behind the AM46 machinery, currently on the order of $30{+}50$ CPU-years. An accepted $m=4$ cut at $n=45$ shrinks the pointed-gluing lists directly through any of the three frozen routes: eliminating degree-20 vertices (a whole gluing side), eliminating the deficiency-$8$ catalog layer (the deepest enumeration tier), or cutting the global deficiency budget $360\to 315$. A decisive negative is still publishable-grade methodology evidence: an exact $q$-refined LP cone with exact rational affine certificates and a disjoint checker is not in print (closest precedent: MR97 §5, $\mathrm{LP}(s,t,n)$, hand-derived and uncertified), and the $m=4$ identity has never before sat inside an exact LP certificate at $n\le 48$ — it was used once, by hand, at $n=49$ in 1997. This plan therefore does **not** promise a cut; it promises an exact, checkable decision.

Known slack from the completed $m=3$ campaign (`data/higher_identity_m3.json`, disposition `M3_NO_CUT_IN_FROZEN_CONE`, 3,215 states, all three route certificates and witnesses exact tight equalities):

| Route | $m{=}3$ exact bound | Acceptance edge |
|---|---|---|
| `total_deficiency` | $360$ | $\le 315$ |
| `degree20_count` | $14310/349\approx 41.00$ | $<1$ |
| `deficiency_ge8_count` | $45$ | $<1$ |

The $m=4$ row adds exactly **one** global equality ($\sum_v h_v=0$, the $m{=}4$ vertex row summed) and therefore two affine certificate coefficients per route ($\varepsilon$ on $h_{lo}$, $\zeta$ on $h_{hi}$). The measured prior comes from `notes/cone_lp_robustness_2026-08-17.md` §4: modeling the added equality as a g-window shrink moves route 2 only ($41.00 \to 36.97$ halved $\to 34.59$ exact-midpoint; exact g-sign $\to 30.79$) and leaves routes 1 and 3 at **zero** reduction even under exact signs. That proxy cannot see the genuinely new per-state interval fields ($t$, diamond, $K_3{+}K_1$, $i_4(\overline{D})$), which is the entire justification for executing this plan — but no route-flip is promised.

Fallback order, recorded on a negative disposition and executed only under a separately approved design: (1) VeriPB certificates for the existing gluing runs; (2) Engström mixed-quadratic tier; (3) complementary SAT on $\mathrm{srg}(45,22,10,11)$-type configurations.

**Architecture:** First implement and exhaustively validate the universal $m=4$ identity kernel, including the MR97 $n=49$ ground-truth replay on the frozen order-24 census. Then stream every published-complete $R(4,5)$ catalog at orders 17–24 once, recording exact per-stratum histograms of $q$, $t$, diamond, $s(K_3{+}K_1)$, and $i_4$ (the complement-side $i_4$ of $(4,5)$-graphs is the new outer-window content), and derive sound outer windows for missing classes from the verified universal 4-vertex relation system. Extend the frozen $n=45$ state cone with $m=4$ row intervals and search for exact affine certificates over the enlarged basis; every certificate carries both the $m=3$ and $m=4$ rows. A separate checker, disjoint from all producers, reconstructs everything from the catalogs.

**Tech Stack:** Python 3.14.3, standard-library `unittest`, integer bitsets, `fractions.Fraction`, SciPy 1.18.0 `linprog` for coefficient discovery only, validated graph6 catalogs, existing frozen modules `check_ramsey`, `subgraph_identities`, `m3_deficiency_cone`.

## Trust Roots

The four $m{=}3$ trust roots carry over verbatim (statements recorded in `data/higher_identity_m3.json:trust_roots`):

1. Each selected order-17..23 file of `data/r45extreme` is the complete census of one fixed-edge $R(4,5)$ class (Angeltveit–McKay census, arXiv:2409.15709); not re-proved locally.
2. The 13 catalogs `r35_1.g6..r35_13.g6` are the complete censuses of $R(3,5)$ graphs (McKay data page); not re-proved locally.
3. `r45_24.g6` is the complete census of the 352,366 $R(4,5,24)$ graphs (arXiv:1703.08768); not re-proved locally.
4. `r55_42some.g6` holds 328 published $R(5,5,42)$ representatives; replay proves the identity on them, not that they exhaust $R(5,5,42)$.

New $m{=}4$ trust roots:

5. **The $m=4$ identity itself**: locally re-derived (module docstring counting ledger, Task 1), exhaustively verified on all labeled graphs $n\le 6$, spot-verified $7$–$8$, replayed on all 656 $n{=}42$ graphs and complements, and replayed against the MR97 Thm 3.2 $n=49$ numbers. Primary-source cross-check anchors: MR97 Thm 2.2 ($m=4$ case) and Thm 3.2 ($n=49$), https://users.cecs.anu.edu.au/~bdm/papers/r55.pdf.
6. **The local motif-window tables** in `data/higher_identity_m4.json`: derived only from hash-pinned frozen catalogs plus the relation system R1–R5 below (whose constants are unit-tested); consumed downstream only through the artifact.

## Global Constraints

- Canonical design: `math/r55/bench/spec.json:next_research_campaign`; this plan is its $m=4$ slice and must not start before the $m{=}3$ checker task has landed (it owns the spec slot and `check_m3_certificate.py` pattern this plan mirrors).
- AI and SciPy remain outside the trust base; every accepted coefficient is rationalized and checked exactly.
- Reuse `check_ramsey` for graph6 parsing and `popcount`, and `subgraph_identities` for `induced_subgraph`, `complement`, `edge_count`, `triangle_count`; do not create a second parser or second triangle counter.
- The frozen $m{=}3$ artifact `data/higher_identity_m3.json` is read-only evidence; $m=4$ results go to a fresh artifact.
- Corrected-math rule: the exact $i_4$ and diamond counters specified in Task 1 supersede the research notes' shortcut forms. The shortcut $i_4\approx\frac16\sum_{uv\notin E}\binom{c_{uv}-2}{2}$ is **provably inexact** (it equals $i_4+X/6$, $X$ = induced one-edge quartets; mismatches on 184/300 random graphs). Fast-versus-naive cross-validation is a hard test requirement.
- Missing inputs, count drift, identity mismatch, or certificate mismatch must exit nonzero.
- Do not change the frozen benchmark engine or candidate-3 artifacts; do not add dependencies.
- No $m=5$, Engström mixed-quadratic, SDP, gluing-engine, or VeriPB implementation work in this plan; the fallback order is recorded, not executed.
- No commit, push, or external mutation without separate authorization; this project tree is currently untracked by Git.
- A valid negative result is `M4_NO_CUT_IN_FROZEN_BASIS`, backed by exact primal witnesses for every instantiable route; never relabel a numerical trend, unresolved certification, or relaxation result as a theorem about all possible uses of $m=4$.
- Do not claim $R(5,5)\le$ anything new.

---

### Task 1: Exact $m=4$ identity kernel

**Files:**
- Create: `math/r55/src/m4_subgraph_identities.py`
- Create: `math/r55/tests/test_m4_subgraph_identities.py`

**Interfaces:**
- Consumes: `check_ramsey.popcount`, `subgraph_identities.{induced_subgraph, complement, edge_count, triangle_count}`.
- Produces:
  - `diamond_count(adj: list[int]) -> int` — induced $K_4$ minus an edge
  - `induced_k3k1_count(adj: list[int]) -> int` — induced $K_3\sqcup K_1$
  - `independent_quad_count(adj: list[int]) -> int` ($=i_4$)
  - `k4_count(adj: list[int]) -> int`
  - `motif4(adj: list[int]) -> Motif4`
  - `m4_vertex_row_scaled(n: int, neighborhood: list[int], dual: list[int]) -> int`
  - `m4_residual_specialized(adj: list[int]) -> int`
  - `m4_residual_universal(adj: list[int]) -> int`

Two identities are stored, both after multiplication by twelve.

Universal (every graph; this is the mechanically exhaustive object):

$$
\sum_v\Big[3(n{+}2{-}2d_v)\,t(N_v)+4\,\mathrm{diamond}(N_v)-6\,s(K_3{+}K_1,N_v)\Big]
+12\sum_v K_4(N_v)-12\sum_v K_4(D_v)=0,
$$

where $K_4(D_v)=i_4(\overline{D_v})$. Specialized to graphs whose every neighborhood is $K_4$-free (in particular every $(5,5)$-graph), the $K_4(N_v)$ term vanishes and the research closed form is recovered with $h_v := 3(n{+}2{-}2d_v)t(N_v)+4\,\mathrm{diamond}(N_v)-6\,s(K_3{+}K_1,N_v)-12\,i_4(\overline{D_v})$ and $\sum_v h_v=0$. The $n=45$ usable row uses $3(47-2d)\in\{21,15,9,3,-3\}$ for $d\in\{20,\dots,24\}$ — note the sign flip at $d=24$.

The module docstring must include the independent counting derivation, in the ledger form of these verified anchor lemmas (all verified in-session during plan drafting; Task 1 re-proves each by test):

- $12\sum_v K_4(D_v)=12\cdot\#\{\text{induced }K_4\sqcup K_1\}$ — each pattern counted once by its isolated vertex;
- $\sum_v t(N_v)=4\,k_4(G)$ (exhaustively verified $n\le 5$; a paw contributes nothing because the pendant cannot contain the triangle);
- $\sum_v \mathrm{diamond}(N_v)=\sum_{S\cong\text{diamond}}|\mathrm{cn}(S)|$, and identically for $s(K_3{+}K_1,N_v)$ and $K_4(N_v)$, where $\mathrm{cn}(S)$ is the set of common neighbors of the 4-set $S$;
- writing $3(n{+}2{-}2d_v)=3(n{-}4)+6(3{-}d_v)$ and regrouping $\sum_v d_v\,t(N_v)=\sum_{S\cong K_4}\sum_{v\in S}d_v$, every remaining term is an incidence count between 4-sets and outside vertices; the ledger closes pattern-by-pattern, with the $K_5$-type incidence carried exactly by the $+12\sum_v K_4(N_v)$ correction.

The exhaustive labeled-graph test is the second, mechanically independent check of this derivation.

Exact fast counters (production path; each cross-validated against the naive 4-subset enumeration in tests):

$$
\mathrm{diamond}(G)=\sum_{uv\notin E} e\big(G[A_u\cap A_v]\big),
\qquad
i_4(G)=\tfrac16\sum_{uv\notin E}\Big[\tbinom{c_{uv}}{2}-e\big(G[W_{uv}]\big)\Big],
$$

where $W_{uv}=\overline{A_u}\cap\overline{A_v}\setminus\{u,v\}$ is the common non-neighborhood, $c_{uv}=|W_{uv}|$, and $s(K_3{+}K_1,G)=\sum_{\text{triangles }T}|\overline{A_a}\cap\overline{A_b}\cap\overline{A_c}\setminus T|$ per triangle $T=\{a,b,c\}$. The diamond formula counts each diamond once via its unique non-edge (adjacent common-neighbor pair); the $i_4$ formula subtracts the induced edges inside $W_{uv}$ — this subtraction is what the research shortcut omits.

Ground-truth replay fixture (MR97 Thm 3.2, $n=49$; all numbers reproduced on `data/r45_24.g6` during plan drafting): exactly two 132-edge $R(4,5,24)$ graphs $X_1,X_2$; both have $t=176$, $\mathrm{diamond}=792$, $s(K_3{+}K_1)=528$; $i_4\in\{138,144\}$; the $n=49,d=24$ row is $9\cdot176+4\cdot792-6\cdot528=1584=12\cdot 132$, forcing per-vertex mean $i_4=132$, while both admissible neighborhood types have $i_4\ge 138>132$ — the identity alone refutes every 24-regular hypothetical Ramsey$(5,5,49)$ configuration.

- [x] **Step 1: Write the exhaustive failing tests**

Create `test_m4_subgraph_identities.py` with a local unlabeled-graph constructor, a naive reference motif counter (4-subset enumeration, test-only), and these observable contracts:

```python
import itertools
import pathlib
import random
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from check_ramsey import parse_graph6_line
from subgraph_identities import complement, induced_subgraph
from m4_subgraph_identities import (  # expected to fail before Task 1 implementation
    diamond_count, independent_quad_count, induced_k3k1_count, k4_count,
    m4_residual_specialized, m4_residual_universal, motif4,
)


def graph_from_mask(n, mask):
    adj = [0] * n
    for bit, (u, v) in enumerate(itertools.combinations(range(n), 2)):
        if (mask >> bit) & 1:
            adj[u] |= 1 << v
            adj[v] |= 1 << u
    return adj


class M4IdentityTests(unittest.TestCase):
    def test_fast_counters_match_naive_through_six_vertices(self):
        for n in range(1, 7):
            for mask in range(1 << (n * (n - 1) // 2)):
                adj = graph_from_mask(n, mask)
                self.assertEqual(diamond_count(adj), naive_diamond(adj))
                self.assertEqual(induced_k3k1_count(adj), naive_k3k1(adj))
                self.assertEqual(independent_quad_count(adj), naive_i4(adj))
                self.assertEqual(k4_count(adj), naive_k4(adj))

    def test_m4_universal_identity_exhaustive_through_six_vertices(self):
        for n in range(1, 7):
            for mask in range(1 << (n * (n - 1) // 2)):
                self.assertEqual(m4_residual_universal(graph_from_mask(n, mask)), 0)

    def test_specialized_identity_on_k4free_neighborhood_graphs(self):
        tested = 0  # plan-drafting reference count over all labeled n<=6: 33,694
        for n in range(1, 7):
            for mask in range(1 << (n * (n - 1) // 2)):
                adj = graph_from_mask(n, mask)
                if any(k4_count(induced_subgraph(
                        adj, [u for u in range(n) if (adj[v] >> u) & 1]))
                       for v in range(n)):
                    continue
                self.assertEqual(m4_residual_specialized(adj), 0)
                tested += 1
        self.assertGreater(tested, 30000)

    def test_specialized_raises_on_k4_neighborhood(self):
        k5 = graph_from_mask(5, (1 << 10) - 1)
        with self.assertRaises(ValueError):
            m4_residual_specialized(k5)

    def test_mutated_coefficients_are_detected(self):
        # Each perturbation is caught separately: diamond coefficient 4 -> 3
        # and K4-correction 12 -> 11 must each yield a nonzero residual on at
        # least one labeled graph on <= 6 vertices; pin explicit witness
        # graphs found at implementation time in place of this scan.
        for mutation in ("diamond_coef", "k4_correction"):
            worst = max(
                abs(mutated_universal_residual(
                    graph_from_mask(n, mask), mutation))
                for n in range(4, 7)
                for mask in range(1 << (n * (n - 1) // 2))
            )
            self.assertGreater(worst, 0, mutation)

    def test_mr97_n49_ground_truth_replay(self):
        graphs = []
        with (ROOT / "data" / "r45_24.g6").open(encoding="ascii") as fh:
            for line in fh:
                if not line.strip():
                    continue
                _, adj = parse_graph6_line(line)
                if sum(row.bit_count() for row in adj) // 2 == 132:
                    graphs.append(adj)
        self.assertEqual(len(graphs), 2)
        for x in graphs:
            f = motif4(x)
            self.assertEqual(f.triangles, 176)
            self.assertEqual(f.diamonds, 792)
            self.assertEqual(f.k3k1, 528)
        self.assertEqual(sorted(independent_quad_count(x) for x in graphs),
                         [138, 144])
        for x in graphs:
            row = 3 * (49 + 2 - 2 * 24) * 176 + 4 * 792 - 6 * 528
            self.assertEqual(row, 1584)          # = 12 * 132, forced i4 mean
            self.assertGreater(independent_quad_count(x), 132)

    def test_all_656_known_r55_graphs_have_zero_specialized_residual(self):
        published = []
        with (ROOT / "data" / "r55_42some.g6").open(encoding="ascii") as fh:
            for line in fh:
                if line.strip():
                    _, adj = parse_graph6_line(line)
                    published.append(adj)
        self.assertEqual(len(published), 328)
        tested = 0
        for adj in published:
            self.assertEqual(m4_residual_specialized(adj), 0)
            self.assertEqual(m4_residual_specialized(complement(adj)), 0)
            tested += 2
        self.assertEqual(tested, 656)

    def test_spot_random_seven_eight_vertices(self):
        rng = random.Random(20260817)
        for _ in range(40):
            n = rng.choice([7, 8])
            adj = [0] * n
            for u, v in itertools.combinations(range(n), 2):
                if rng.random() < rng.choice([0.3, 0.5, 0.7]):
                    adj[u] |= 1 << v
                    adj[v] |= 1 << u
            self.assertEqual(m4_residual_universal(adj), 0)
```

`mutated_universal_residual(adj, mutation)` (mutation selector `"diamond_coef"` perturbs the diamond coefficient 4->3; `"k4_correction"` perturbs the correction 12->11) and the `naive_*` counters live in the test file (the naive ones as 4-subset enumerations; the mutated one recomputing the universal residual with the stated coefficient perturbed).

- [x] **Step 2: Run the test and observe the intended failure**

Run:

```bash
cd /Users/jinleic/jinleic-workspace/math
./.venv/bin/python -m unittest r55/tests/test_m4_subgraph_identities.py -v
```

Expected: import failure for `m4_subgraph_identities`.

- [x] **Step 3: Implement the exact kernel**

Use this implementation shape in `m4_subgraph_identities.py`:

```python
#!/usr/bin/env python3
from collections import namedtuple

from check_ramsey import popcount
from subgraph_identities import complement, induced_subgraph, triangle_count

Motif4 = namedtuple(
    "Motif4", "order triangles diamonds k3k1 independent_quads cliques"
)


def _pairs_inside(adj, mask):
    """Induced edge count inside the vertex bitmask, via per-vertex rows."""
    edges = 0
    while mask:
        bit = mask & -mask
        v = bit.bit_length() - 1
        edges += popcount(mask & adj[v])
        mask ^= bit
    return edges // 2


def diamond_count(adj):
    # each induced diamond: unique non-edge (its degree-2 pair) whose
    # common-neighbor set contains exactly the adjacent pair {w, x}
    n = len(adj)
    total = 0
    for u in range(n):
        for v in range(u + 1, n):
            if (adj[u] >> v) & 1:
                continue
            total += _pairs_inside(adj, adj[u] & adj[v])
    return total


def independent_quad_count(adj):
    # (1/6) * sum over non-edges uv of [C(c,2) - e(W)],
    # W = common non-neighbors of u and v, u and v excluded; the edge
    # subtraction is REQUIRED: without it the sum equals 6*i4 + X
    # where X counts induced one-edge quartets (each contributes 1).
    n = len(adj)
    dual = complement(adj)
    total = 0
    for u in range(n):
        for v in range(u + 1, n):
            if (adj[u] >> v) & 1:
                continue
            w = dual[u] & dual[v] & ~((1 << u) | (1 << v))
            c = popcount(w)
            total += c * (c - 1) // 2 - _pairs_inside(adj, w)
    assert total % 6 == 0
    return total // 6


def k4_count(adj):
    return independent_quad_count(complement(adj))


def induced_k3k1_count(adj):
    n = len(adj)
    dual = complement(adj)
    total = 0
    for a in range(n):
        later = adj[a] & ~((1 << (a + 1)) - 1)
        while later:
            b1 = later & -later
            b = b1.bit_length() - 1
            tri = adj[a] & adj[b] & ~((1 << (b + 1)) - 1)
            while tri:
                b2 = tri & -tri
                c = b2.bit_length() - 1
                total += popcount(
                    dual[a] & dual[b] & dual[c]
                    & ~((1 << a) | (1 << b) | (1 << c))
                )
                tri ^= b2
            later ^= b1
    return total


def _neighborhood(adj, v):
    return induced_subgraph(
        adj, [u for u in range(len(adj)) if (adj[v] >> u) & 1]
    )


def _dual(adj, v):
    n = len(adj)
    return induced_subgraph(
        adj, [u for u in range(n) if u != v and not ((adj[v] >> u) & 1)]
    )


def m4_vertex_row_scaled(n, neighborhood, dual):
    if len(neighborhood) + len(dual) != n - 1:
        raise ValueError("neighborhood and dual orders must sum to n-1")
    d = len(neighborhood)
    return (
        3 * (n + 2 - 2 * d) * triangle_count(neighborhood)
        + 4 * diamond_count(neighborhood)
        - 6 * induced_k3k1_count(neighborhood)
        - 12 * independent_quad_count(complement(dual))
    )


def m4_residual_specialized(adj):
    n = len(adj)
    for v in range(n):
        if k4_count(_neighborhood(adj, v)):
            raise ValueError("specialized identity needs K4-free neighborhoods")
    return sum(
        m4_vertex_row_scaled(n, _neighborhood(adj, v), _dual(adj, v))
        for v in range(n)
    )


def m4_residual_universal(adj):
    n = len(adj)
    return sum(
        m4_vertex_row_scaled(n, _neighborhood(adj, v), _dual(adj, v))
        + 12 * k4_count(_neighborhood(adj, v))
        - 12 * k4_count(_dual(adj, v))
        + 12 * independent_quad_count(complement(_dual(adj, v)))
        for v in range(n)
    )
```

Note the universal residual simplifies to
$\sum_v[\text{row}(N_v)+12K_4(N_v)-12K_4(D_v)]$ since $K_4(D_v)=i_4(\overline{D_v})$; the displayed shape keeps the row function single-sourced. Keep whatever explicit form makes the docstring ledger match the code one-to-one.

- [x] **Step 4: Run exhaustive tests and the existing graph-parser selftest**

Run:

```bash
cd /Users/jinleic/jinleic-workspace/math
./.venv/bin/python -m unittest r55/tests/test_m4_subgraph_identities.py -v
./.venv/bin/python r55/src/check_ramsey.py --selftest
./.venv/bin/python -m unittest r55/tests/test_subgraph_identities.py -v
```

Expected: all new tests `OK` (including the 656-graph replay and the $n=49$ replay); existing suites unchanged and `selftest OK`.

---

### Task 2: Exact motif windows and the extended $n=45$ state cone

**Files:**
- Create: `math/r55/src/m4_deficiency_cone.py`
- Create: `math/r55/tests/test_m4_deficiency_cone.py`
- Create at execution time: `math/r55/data/higher_identity_m4.json`

**Interfaces:**
- Consumes: Task 1 counters, `m3_deficiency_cone.{load_edge_bounds, r35_edge_windows, degree_histograms, q_outer_interval, q_interval, budget2, degree_bounds}` (frozen, verified $m{=}3$ infrastructure — one owner per fact), `data/structural_tables.json`, `data/VALIDATION.json`, `data/VALIDATION_extreme.json`, the 36 validated edge-extremal `r45extreme/*.g6` catalogs at orders 17–23, and the complete `r45_24.g6` census.
- Produces:
  - `motif4_vector(adj) -> tuple[int, int, int, int, int]` ($t$, $q$, diamond, $s(K_3{+}K_1)$, $i_4$; one streaming pass)
  - `stream_catalog_motifs(path) -> dict[tuple[int, int], dict[str, dict[int, int]]]`
  - `catalog_motif_windows(records) -> dict[tuple[int, int], dict[str, tuple[int, int]]]`
  - `outer_motif_window(order, edges, windows, edge_bounds) -> dict[str, tuple[int, int]]`
  - `affine_interval(coef, lo, hi) -> tuple[int, int]` (sign-aware: the $t$-coefficient flips at $d=24$)
  - `m4_state_interval(d, a, b, windows, edge_bounds, catalog_windows) -> State`
  - `build_states(windows, edge_bounds, catalog_windows) -> list[State]`
  - deterministic analysis JSON with source hashes, both replays, exact per-catalog motif histograms, sound missing-class windows, and every state interval.

Per-graph budget (the whole-task cost driver). For one $R(4,5)$ graph of order $r$ with $e$ edges and $N=\binom r2-e$ non-edges:

| Piece | Operations |
|---|---|
| degrees, $q$ | $r$ popcounts $+$ triangle pass ($e$ AND+popcount) |
| diamond | per non-edge: 1 AND + popcount + $|A_u\cap A_v|$ row-ANDs $\Rightarrow N\cdot(2+\bar c)$ |
| $i_4$ | per non-edge: 2 ANDs + popcount + $|W_{uv}|$ row-ANDs $\Rightarrow N\cdot(3+\bar w)$ |
| $s(K_3{+}K_1)$ | per triangle: 2 ANDs + 1 popcount $\Rightarrow 3t$ |

Measured anchors on this machine (Apple M3 Ultra, pure Python, plan-drafting run): $0.287$ ms/graph (order 19, $e=90$), $0.421$ ms/graph (order 23, $e=119$), $0.507$ ms/graph (order 24). Catalog total $8{,}500{,}211$ graphs $= 8{,}147{,}845$ in the 36 validated order-17..23 files $+ 352{,}366$ in `r45_24.g6`, matching the $m{=}3$ artifact input records exactly. Projected single-process wall time $\approx 0.8$ h; the fail-closed budget is **2 h single-process**. Streaming discipline: read one graph6 line at a time, fold into per-`(order, edges)` histograms, retain no graph, hash every file first, and never hold more than one adjacency list. Multiprocessing is explicitly unnecessary; do not add it.

For a missing class $(r,e)$ the sound windows for diamond, $s(K_3{+}K_1)$, $i_4$ come from the universal 4-vertex relation system over the 11 induced motifs (order $E_4,K_2{+}2K_1,2K_2,P_3{+}K_1,P_4,K_{1,3},K_3{+}K_1,C_4,\text{paw},\text{diamond},K_4$):

$$
\begin{aligned}
R1&:\ \textstyle\sum_J s(J)=\binom r4 &
R2&:\ \textstyle\sum_J e(J)\,s(J)=e\binom{r-2}2\\
R3&:\ \textstyle\sum_J e(J)^2 s(J)=e\binom{r-2}2+2W(r-4)+e(e-1) &
R4&:\ \textstyle\sum_J \big(\sum_{v\in J}\deg_J(v)^2\big)s(J)=\sum_x\Big[d_x\binom{r-1-d_x}{2}+4\binom{d_x}2(r-1-d_x)+9\binom{d_x}3\Big]\\
R5&:\ \textstyle\sum_J \tau(J)\,s(J)=(r-3)\,t &
\end{aligned}
$$

with $W=\sum_x\binom{d_x}2$ and $\tau$ the triangle multiplicity ($K_3{+}K_1$:1, paw:1, diamond:2, $K_4$:4). $R1$–$R5$ verified on 60 random graphs with exactly these constants. Three corrections to the research notes are baked in: $R3$'s $p_3+3t$ collapses to $W$; $R4$'s right side is the degree-histogram form (the quoted $2e\binom{r-2}2$ is false); $R5$ needs multiplicities. The system has rank 5 (6 free ordinals), and **no single motif is point-determined**: $s$-vectors with identical aggregates differ in every coordinate (verified: nullspace touches all 11; $(r,e,p_3,t)$-twins vary in all 11). Windows therefore come from an exact optimization, not from solving for a coordinate.

**RHS intervals (soundness-critical).** Each relaxed right-hand side is the *exact min and max of its aggregate over the same Ramsey-feasible integer degree histograms already enumerated by the imported `q_outer_interval` machinery*: R2's RHS $e\binom{r-2}2$ is class-determined and stays **exact**; R3's aggregate is $W$; R4's aggregate is $\sum_x f(d_x)$ with $f(d)=d\binom{r-1-d}2+4\binom d2(r-1-d)+9\binom d3$; R5's aggregate is $t\in[W-q_{hi},\,W-q_{lo}]$ from the $m{=}3$ $q$-machinery. No coarser proxy is permitted: an RHS interval that excludes a feasible histogram's true value is the unsound direction (windows would come out too narrow and a certificate checked against them could fail on real graphs).

**Exact acceptance rule.** The window of each needed motif over $\{s\in\mathbb{R}_{\ge0}^{11}:\ R1\ \text{exact},\ R2\ \text{exact},\ R3\text{–}R5\ \text{two-sided}\}$ is computed by exhaustive enumeration of basic feasible solutions in `Fraction` arithmetic (holding R1 *and* R2 exact: 11 nonnegativity + 6 equality + 6 two-sided rows keeps candidate bases at $\binom{17}{9}\approx 2.4\times10^4$ per class). SciPy output is a **pruning hint only**: an endpoint is accepted only when it is the exact min/max over *all* `Fraction`-verified basic solutions. Verifying a SciPy-proposed vertex alone never certifies optimality. On any certification failure the class falls back to the trivially sound $[0,\binom r4]$.

**Envelope-phase budget.** The Step 1 tests pin one fixed class's deterministic verified-basic-solution count and record its measured wall anchor (informational, this machine). Fail-closed ceiling for the whole envelope phase — 92 missing strata plus 53 audit classes — is **30 minutes single-process**; exceeding it is a hard nonzero exit, never a silent widening. Present classes keep their exact catalog windows; the LP windows are additionally computed for present classes as a self-audit (see Open Questions).

Each state $(d,a,b)$ extends the frozen $m{=}3$ record with the $m{=}4$ row interval, where $x$ is the $R(4,5,d)$ neighborhood with $e(x)=E_{45}(d)-a$ and $y=\overline{D_v}$ the $R(4,5,m)$ graph with $e(y)=E_{45}(m)-b$, $m=44-d$:

$$
h \in \Big[\,\mathrm{row}_{lo},\ \mathrm{row}_{hi}\,\Big],\qquad
\mathrm{row}=3(47-2d)\,t(x)+4\,\mathrm{diamond}(x)-6\,s(K_3{+}K_1,x)-12\,i_4(y),
$$

computed sign-aware via `affine_interval` over the four stratum windows, and $\sum_v h_v=0$ gives $\sum_v h_{lo}\le 0\le\sum_v h_{hi}$ — the one new global pair of inequalities the $m{=}4$ row contributes.

- [x] **Step 1: Write failing envelope tests**

Tests must cover (naive reference counters allowed as test helpers):

1. fast `motif4_vector` equals the naive 4-subset vector on 250 consecutive graphs of `r45extreme/r4517.77.g6` and on the two 132-edge order-24 graphs;
2. `affine_interval` sign handling: `affine_interval(-3, (10, 20)) == (-60, -30)`;
3. every relation $R1$–$R5$ holds with the exact constants above on 60 seeded random graphs ($n\le 7$);
4. exact catalog windows override outer windows, for each of the five motifs;
5. for the three fixture stratum pairs of the $m{=}3$ tests — $(20,100){\times}(24,132)$, $(21,107){\times}(23,122)$, $(22,114){\times}(22,114)$ — the exact $h$-value of every catalog-graph pair (first 20 of each side) lies inside `m4_state_interval(...).h_lo, h_hi`, and likewise the $m{=}3$ $g$-interval still contains the $m{=}3$ expression;
6. `outer_motif_window` returns windows containing every catalog graph's motif value for at least three present classes (audit direction);
7. the state count matches the $m{=}3$ enumeration domain ($d\in\{20,\dots,24\}$, full $(a,b)$ boxes; see Open Questions) and every state satisfies `h_lo <= h_hi`;
8. the 656-graph specialized replay and the $n=49$ replay both re-run through the Task 1 kernel.

- [x] **Step 2: Run tests and observe missing-interface failures**

Run:

```bash
cd /Users/jinleic/jinleic-workspace/math
./.venv/bin/python -m unittest r55/tests/test_m4_deficiency_cone.py -v
```

Expected: import failures for the cone interfaces.

- [x] **Step 3: Implement streaming windows and the extended State**

Use exact integers; extend (do not mutate) the frozen $m{=}3$ shape:

```python
@dataclass(frozen=True, order=True)
class State:
    d: int
    a: int
    b: int
    deficiency: int
    excess_balance: int
    g_lo: int          # m=3 row, reused unchanged from the m=3 machinery
    g_hi: int
    h_lo: int          # m=4 row, new
    h_hi: int
```

`build_states` reuses `m3_deficiency_cone.m3_state_interval` logic for the $g$ endpoints (import, do not copy) and adds the $h$ endpoints from the four motif windows of the two strata. It may omit a pair only when outer degree-histogram constraints prove one side infeasible, exactly as in Task 2 of the $m{=}3$ plan; missing $R(4,5)$ catalogs are never treated as infeasibility.

- [x] **Step 4: Add fail-closed input provenance and the v2 artifact**

The CLI must:

1. hash every input file with SHA-256;
2. load edge windows from the validated complete `r35_k.g6` records (reuse `r35_edge_windows`);
3. assert all 328 lines plus complements give 656 zero specialized residuals;
4. replay the $n=49$ ground truth (two 132-edge graphs, row $1584=12\cdot132$, $i_4$ multiset $\{138,144\}$);
5. stream every validated order-17..23 file and the complete order-24 census, one line at a time, recording exact per-class histograms of $q$, $t$, diamond, $s(K_3{+}K_1)$, $i_4$;
6. verify every observed motif value lies inside its `outer_motif_window` (envelope self-audit);
7. use exact catalog windows for present classes, LP windows for missing classes, trivial windows only on LP certification failure;
8. write JSON only after all checks pass, sorted keys, deterministic lists, atomic replace via temporary sibling.

The JSON schema is fixed before search:

```text
schema_version: 2
campaign_id: higher_order_identity_positive_deficiency_m4
inputs: sorted records {relative_path, sha256, bytes, graph_count}
trust_roots: the six statements of this plan's Trust Roots section
replay: {published_graphs: 328, complements: 328, residual_zero: 656}
n49_replay: {extremal_graphs: 2, row_scaled: 1584, forced_i4_mean: 132,
             observed_i4: [138, 144]}
r35_edge_windows: sorted records {order, edge_min, edge_max}
catalog_motif_histograms: sorted records
  {order, edges, source, graph_count,
   histograms: {q: [{value, count}], t: [...], diamond: [...],
                k3k1: [...], i4: [...]}}
catalog_motif_windows: sorted records
  {order, edges, source, windows:
   {q: [lo, hi], t: [...], diamond: [...], k3k1: [...], i4: [...]}}
outer_motif_windows: sorted records {order, edges, method, windows}
states: sorted records
  {d, a, b, deficiency, excess_balance, g_lo, g_hi, h_lo, h_hi,
   x_interval_source, y_interval_source,
   x_motif_source, y_motif_source}
searches: [] until Task 3
disposition: M4_ANALYZED until Task 3
```

No wall-clock timestamp or absolute path enters the canonical JSON.

- [x] **Step 5: Run focused tests and full analysis**

Run:

```bash
cd /Users/jinleic/jinleic-workspace/math
./.venv/bin/python -m unittest r55/tests/test_m4_deficiency_cone.py -v
./.venv/bin/python r55/src/m4_deficiency_cone.py \
  --output r55/data/higher_identity_m4.json
```

Expected terminal contract:

```text
M4 REPLAY: 656/656 specialized residual zero
M4 N49 REPLAY: extremal=2 row=1584 forced_i4_mean=132 observed=[138,144]
M4 STREAM: 8500211 graphs in <seconds> (budget 7200)
M4 ENVELOPES: catalog violations=0 outer-window audits=ALL
M4 STATES: N states (N >= 1)
WROTE r55/data/higher_identity_m4.json
```

---

### Task 3: Exact affine certificates with $m=3$ and $m=4$ rows

**Files:**
- Create: `math/r55/src/search_m4_cuts.py`
- Create: `math/r55/tests/test_search_m4_cuts.py`
- Modify at execution time: `math/r55/data/higher_identity_m4.json`

**Interfaces:**
- Consumes: Task 2 states and the frozen objective registry.
- Produces an immutable `Certificate` with fields `objective_id`, `alpha`, `beta`, `gamma`, `delta`, `epsilon`, `zeta`; rational values serialize as reduced fraction strings.
- Produces an immutable `PrimalWitness` with `objective_id` and sorted `(state_index, weight)` rational pairs.
- Produces `exact_affine_bound(states, objective_id) -> Certificate`, `verify_certificate(states, certificate)`, `verify_primal_witness(states, witness)`, and `verify_search_record(states, record)`.
- Every certificate satisfies $\gamma,\delta,\varepsilon,\zeta\ge0$ and, for every state,

$$
objective(state)\le \alpha+\beta\,balance(state)+\gamma\, g_{lo}(state)-\delta\, g_{hi}(state)
 +\varepsilon\, h_{lo}(state)-\zeta\, h_{hi}(state).
$$

The cone basis is the frozen $m{=}3$ basis plus the $m{=}4$ row: every state carries both intervals, every certificate format carries all six coefficients, and the dual discovery optimizes over all of them ($\gamma$ or $\delta$ being zero is an output, never a constraint — the $m{=}3$ rows are required and reused). Summing over 45 vertices uses $\sum balance=0$, $\sum g_{lo}\le0\le\sum g_{hi}$, and the new $\sum h_{lo}\le0\le\sum h_{hi}$, proving the exact global bound $\sum objective\le45\alpha$.

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

from m4_deficiency_cone import State
from search_m4_cuts import (
    PrimalWitness, exact_affine_bound, verify_certificate,
    verify_primal_witness,
)


class ExactCutTests(unittest.TestCase):
    def test_rationalized_certificate_checks_every_state(self):
        states = [
            State(20, 0, 0, 0, -12, -4, 2, -8, 1),
            State(24, 0, 0, 0, -12, -2, 4, -6, 3),
        ]
        cert = exact_affine_bound(states, "total_deficiency")
        verify_certificate(states, cert)

    def test_tampered_alpha_is_rejected(self):
        states = [
            State(20, 0, 0, 0, -12, -4, 2, -8, 1),
            State(24, 0, 0, 0, -12, -2, 4, -6, 3),
        ]
        cert = exact_affine_bound(states, "total_deficiency")
        tampered = replace(cert, alpha=str(Fraction(cert.alpha) - 1))
        with self.assertRaises(ValueError):
            verify_certificate(states, tampered)

    def test_negative_epsilon_is_rejected(self):
        states = [State(20, 0, 0, 0, 0, -1, 1, -1, 1)]
        cert = exact_affine_bound(states, "deficiency_ge8_count")
        tampered = replace(cert, epsilon="-1")
        with self.assertRaises(ValueError):
            verify_certificate(states, tampered)

    def test_exact_primal_rejection_witness_with_m4_rows(self):
        states = [State(20, 0, 0, 0, 0, -1, 1, -2, 2)]
        witness = PrimalWitness("degree20_count", ((0, "45"),))
        self.assertEqual(verify_primal_witness(states, witness), Fraction(45))

    def test_primal_witness_violating_m4_row_is_rejected(self):
        states = [State(20, 0, 0, 0, 0, -1, 1, 1, 2)]  # h_lo > 0
        witness = PrimalWitness("degree20_count", ((0, "45"),))
        with self.assertRaises(ValueError):
            verify_primal_witness(states, witness)

    def test_wrong_objective_id_is_rejected(self):
        states = [State(20, 0, 0, 0, 0, -1, 1, -1, 1)]
        witness = PrimalWitness("not_registered", ((0, "45"),))
        with self.assertRaises(ValueError):
            verify_primal_witness(states, witness)
```

- [x] **Step 2: Run tests and observe missing-interface failures**

Run:

```bash
cd /Users/jinleic/jinleic-workspace/math
./.venv/bin/python -m unittest r55/tests/test_search_m4_cuts.py -v
```

- [x] **Step 3: Implement exact upper certificates and lower witnesses**

Freeze this registry in source before the first full run — identical routes and edges to $m{=}3$:

```text
total_deficiency: objective=state.deficiency, accept upper_bound <= 315
degree20_count: objective=int(state.d == 20), accept upper_bound < 1
deficiency_ge8_count: objective=int(state.deficiency >= 8), accept upper_bound < 1
```

For an actual 45-vertex graph the two count objectives are nonnegative integers, so an exact upper bound below one proves the count is zero. No integrality assumption is used for the continuous relaxation itself.

Solve the continuous dual discovery LP with variables $(\alpha,\beta,\gamma,\delta,\varepsilon,\zeta)$, bounds `[(None,None), (None,None), (0,None), (0,None), (0,None), (0,None)]`, objective `[1, 0, 0, 0, 0, 0]`, and one row
`[-1, -balance, -g_lo, g_hi, -h_lo, h_hi] <= -objective(state)` per state. Convert the returned $(\beta,\gamma,\delta,\varepsilon,\zeta)$ to `Fraction` using the frozen denominator ladder `(100, 1000, 10000, 100000, 1000000)`, then derive $\alpha$ exactly as the max over states of $objective-\beta\,balance-\gamma\,g_{lo}+\delta\,g_{hi}-\varepsilon\,h_{lo}+\zeta\,h_{hi}$. `verify_certificate` resolves `objective_id` through the frozen registry, rejects unknown IDs, checks the four coefficient signs and every state inequality, and returns $45\alpha$. No caller-provided lambda or stored boolean controls acceptance.

For every route not accepted by an exact upper certificate, solve the primal discovery LP with variables $x_s\ge0$:

```text
maximize sum objective(s) * x_s
subject to sum x_s = 45
           sum excess_balance(s) * x_s = 0
           sum g_lo(s) * x_s <= 0
           sum g_hi(s) * x_s >= 0
           sum h_lo(s) * x_s <= 0
           sum h_hi(s) * x_s >= 0
```

Use HiGHS dual simplex only to identify a sparse basic support (weights above `1e-9`). If exact reconstruction fails, form one deterministic candidate set from those indices and the 24 states with smallest `(absolute_dual_slack, state_index)` — `FALLBACK_SUPPORT = 24`, `MAX_SUPPORT = 6` (2 equalities + 4 inequalities span a basis) — and enumerate subsets of at most six states with every active/inactive choice for the four inequality rows, solving each basis by `Fraction` Gaussian elimination. Do not expand this set after seeing results: failure becomes `CERTIFICATION_UNRESOLVED`.

`verify_primal_witness` independently requires nonnegative weights, total weight 45, zero excess balance, aggregate $g_{lo}\le0\le g_{hi}$, aggregate $h_{lo}\le0\le h_{hi}$, and recomputes the exact objective. A route is `REJECTED_BY_EXACT_WITNESS` only when the witness value is $>315$ for `total_deficiency` or $\ge1$ for either count objective. Otherwise it is `CERTIFICATION_UNRESOLVED`; solver status or a rounded coefficient never proves a negative result.

Route 4 (`required_local_family`, a named $(d,a,b)$ family required by a hash-pinned complete-gluing-cover manifest) stays `UNAVAILABLE_NO_COVER_CERTIFICATE`: no complete $n=45$ gluing-cover manifest exists in this repository. Do not tune objective definitions, thresholds, denominator ladder, or support rules after seeing results.

- [x] **Step 4: Run the cut search and persist exact decisions**

Run:

```bash
cd /Users/jinleic/jinleic-workspace/math
./.venv/bin/python r55/src/search_m4_cuts.py \
  --analysis r55/data/higher_identity_m4.json
```

Expected terminal contract is exactly one of:

```text
M4 CUT SEARCH COMPLETE status=M4_ACCEPTED_CUT accepted=N (N >= 1)
```

```text
M4 CUT SEARCH COMPLETE status=M4_NO_CUT_IN_FROZEN_BASIS accepted=0
```

```text
M4 CUT SEARCH COMPLETE status=CERTIFICATION_UNRESOLVED accepted=0
```

`M4_NO_CUT_IN_FROZEN_BASIS` requires an exact rejection witness for every instantiable route. `CERTIFICATION_UNRESOLVED` exits nonzero and cannot trigger a scientific disposition. Every persisted search record contains its objective ID, exact upper certificate and bound when available, exact lower witness and value when available, derived route status, and no unchecked acceptance field.

---

### Task 4: Disjoint checker and evidence-grade disposition

**Files:**
- Create: `math/r55/src/check_m4_certificate.py`
- Create: `math/r55/tests/test_check_m4_certificate.py`
- Modify: `math/r55/bench/spec.json`
- Modify: `math/PROGRESS.md`
- Modify: `math/r55/README.md`
- Modify: `math/r55/notes/structural_constraints_2026-08-15.md`

This task mirrors the (by then landed) `check_m3_certificate.py` discipline and must not begin before that checker exists.

**Interfaces:**
- Consumes: frozen validation files, graph catalogs, and `data/higher_identity_m4.json`.
- Produces an independent pass/fail verdict. It must not import `subgraph_identities`, `m4_subgraph_identities`, `m3_deficiency_cone`, `m4_deficiency_cone`, or `search_m4_cuts`.

- [x] **Step 1: Write checker corruption tests**

Create a temporary copy of the analysis artifact and independently test rejection of:

1. one changed source SHA-256;
2. one changed state endpoint ($g$ or $h$);
3. one reduced `alpha` coefficient;
4. one changed primal weight;
5. one false exact bound or route status;
6. one changed 656-graph replay count or $n=49$ replay field;
7. one changed motif histogram count (any of the five motifs);
8. one changed outer-window endpoint (any motif).

The untouched artifact must pass. Separately exhaust all labeled graphs through six vertices against the checker's own universal residual implementation, and show that changing the scaled diamond coefficient from 4 to 3, the $i_4$ edge-subtraction, or the $K_4(N_v)$ correction 12 to 11 each fails on an explicit witness graph.

- [x] **Step 2: Implement the independent checker**

The checker deliberately reimplements this small trust kernel:

- graph complement, induced subgraph, edge, triangle, diamond, $K_3{+}K_1$, $K_4$, and $i_4$ counts (exact fast forms, cross-validated against its own naive enumerations);
- the universal and specialized $m{=}4$ residuals, recomputed as zero on all 328 published graphs and their complements, and the $n=49$ replay numbers ($1584=12\cdot132$; $\{138,144\}$);
- degree-histogram enumeration and the $m{=}3$ $q$-window arithmetic (including $q=W-t$ outer intervals) for the $g$ endpoints;
- the relation system $R1$–$R5$ with this plan's corrected constants and the LP outer-window derivation, recomputed for every class the artifact claims;
- exact per-class motif histograms re-streamed from all applicable published-complete files;
- all `(d,a,b)` state endpoints and the five aggregate sign conditions;
- the frozen objective registry, six-coefficient exact dual certificates, exact primal witnesses, bounds, thresholds, and derived route/disposition statuses;
- source hashes, catalog counts, and frozen replay counts.

It may import only Python standard-library modules and `check_ramsey.parse_graph6_line`/`popcount`. It must not trust any producer count, bound, boolean, or disposition without recomputation. Any mismatch prints the first exact counterexample and exits 1.

- [x] **Step 3: Run all focused tests and the independent checker**

Run:

```bash
cd /Users/jinleic/jinleic-workspace/math
./.venv/bin/python -m unittest discover -s r55/tests -p 'test_*m4*.py' -v
./.venv/bin/python r55/src/check_m4_certificate.py \
  r55/data/higher_identity_m4.json
```

Expected: focused tests `OK`; checker prints `M4 EVIDENCE VERIFIED` and the exact disposition.

- [x] **Step 4: Apply the frozen scientific gate**

If at least one exact upper certificate meets a quantified success route:

- set `next_research_campaign.status` to `M4_ACCEPTED_PENDING_NOVELTY_AUDIT`;
- record the theorem statement with exact coefficients and hypotheses in the existing structural note;
- conduct the primary-source novelty audit (MR97 Thm 2.2/3.2 and §5, Engström arXiv:1002.4304, Angeltveit–McKay arXiv:2409.15709) before using "new" or "publishable";
- do not claim $R(5,5)\le45$.

If every instantiable route is rejected by an exact primal witness:

- set status to `M4_NO_CUT_IN_FROZEN_BASIS`;
- record each exact witness and state only that the frozen continuous basis $\{balance, g, h\}$ cannot prove an accepted route;
- preserve all artifacts and hashes;
- record the fallback order (VeriPB certificates $\to$ Engström mixed-quadratic $\to$ $\mathrm{srg}(45,22,10,11)$ complementary SAT) as the decision point for the next separately reviewed campaign; the $T$-family identity line ($m=2,3,4$; $m\ge5$ void) is then closed at $n=45$. The spec's `first_exact_basis` additionally contains the Engström mixed identity, so the basis-level stop condition of `next_research_campaign` remains open until that tier is separately reviewed and executed.

If any route is `CERTIFICATION_UNRESOLVED`:

- set status to `M4_CERTIFICATION_UNRESOLVED`;
- preserve the verified upper/lower bounds without a positive or negative claim;
- stop and record the exact unresolved gap and basis-reconstruction failure.

- [x] **Step 5: Final integrity verification**

Parse the $m{=}4$ artifact, verify every recorded SHA-256, rerun the checker, confirm the frozen $m{=}3$ artifact, benchmark hashes, and candidate artifacts are unchanged, and rerun both $m{=}3$ and $m{=}4$ focused suites. Report the mathematical result first, then files and verification evidence.

---

## Plan Self-Review

- **Spec coverage:** Implements the approved $m=4$ slice of `next_research_campaign`: identity derivation and replay (including the $n=49$ ground truth), exact catalog motif vectors for $q$, $t$, diamond, $K_3{+}K_1$, $i_4$ with sound missing-class envelopes from the corrected relation system, the finite cone extended by the $m=4$ row on top of the reused $m{=}2$/$m{=}3$ rows, exact six-coefficient certificates and witnesses, solver outside the trust base, a disjoint checker, all instantiable quantified routes, and fail-closed positive, negative, or unresolved dispositions with the recorded fallback order.
- **Placeholder scan:** No unresolved implementation placeholders, deferred code bodies, or production stubs remain. The two "pin an explicit witness graph" notes in Task 1 tests are implementation-time constants discovered by the test itself, not stubs.
- **Type consistency:** `State` (nine fields), `Motif4`, `Certificate` (six coefficients), `PrimalWitness`, artifact path `data/higher_identity_m4.json` at `schema_version: 2`, and the CLI terminal contracts are consistent across all four tasks; the $m{=}3$ `State` and artifact are consumed read-only and never mutated.
- **Scope:** No $m=5$, mixed-quadratic, SDP, gluing-engine, VeriPB, benchmark, or raw-generation work enters this plan.

**Novelty table** (from the verified research inputs):

| Candidate deliverable | Closest precedent | Why it would be new |
|---|---|---|
| Exact $q$-refined LP cone with exact rational affine certificates and a disjoint checker | MR97 §5, $\mathrm{LP}(s,t,n)$ — hand-derived, uncertified | machine-certified exact cone over complete published catalogs; not in print |
| $m=4$ identity inside an exact LP certificate at $n\le48$ | MR97 Thm 3.2 — one hand evaluation at $n=49$ (1997) | automated cone plus certificate at $n=45$; never published |
| Decisive negative with exact witnesses | none for this basis | publishable methodology evidence closing the algebraic line |

**If a cut accepts, the theorem is** (route-dependent, exact quantifiers): *For every graph $G$ with $|V(G)|=45$, $\omega(G)\le4$, and $\alpha(G)\le4$, it holds that* — `total_deficiency`: $\sum_{v\in V(G)} s(v)\le 45\alpha\le315$; `degree20_count`: $\min_v d_v\ge21$; `deficiency_ge8_count`: $s(v)\le7$ for every $v$ — *conditional on the named catalog-completeness trust roots, with the recorded rational certificate $(\alpha,\beta,\gamma,\delta,\varepsilon,\zeta)$ verified by the disjoint checker.* None of these claims $R(5,5)\le45$.

**Research-input corrections applied during plan drafting** (all verified computationally in-session; the plan's tests enforce them):

1. $i_4$: the shortcut $\frac16\sum_{uv\notin E}\binom{c_{uv}-2}{2}$ is inexact ($=i_4+X/3$); exact form subtracts induced edges inside the common non-neighborhood.
2. diamond: fast form is $\sum_{uv\notin E} e(G[A_u\cap A_v])$ (adjacent, not all, common-neighbor pairs).
3. $R4$ right side is the degree-histogram form; $R5$ needs triangle multiplicities $(1,1,2,4)$; $R3$'s $p_3+3t$ collapses to $W$.
4. No motif (including $K_3{+}K_1$ and $i_4$) is point-determined by $(r,e,p_3,t)$ or by the five aggregates; missing-class windows must be optimizations over the relation polytope, which Task 2 does.

**Open questions for the main agent:**

- (a) Artifact schema v2: fresh file `data/higher_identity_m4.json` (this plan's choice, keeping the frozen $m{=}3$ artifact intact as evidence) versus extending `higher_identity_m3.json` in place. Recommendation: fresh file.
- (b) State set: re-enumerate the full $(d,a,b)$ domain with the $m{=}3$ feasibility pruning (this plan; expected to reproduce the 3,215-state domain) versus importing the frozen 3,215 records and only adding $h$ endpoints. Recommendation: re-enumerate and assert the count; import-only would silently inherit any future $m{=}3$ changes.
- (c) $i_4$ envelope self-audit scope: compute LP outer windows for all classes (missing authoritative, present audit-only — this plan) versus present classes only. Recommendation: all classes; the audit-only recomputations double as relation-system validation on real data.
