# Engström Mixed Identity $n=45$ Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use `subagent-driven-development` (recommended) or `executing-plans` to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

## Executive summary

**Goal:** Determine, with exact arithmetic and an independent checker, whether the Engström mixed (non-separable) degree-4 neighborhood/anti-neighborhood identity proves a publishable positive-deficiency or degree restriction on every hypothetical Ramsey$(5,5,45)$ graph. A decisive negative (`MIXED_NO_CUT_IN_FROZEN_BASIS`) is an allowed, evidence-grade outcome: it closes the campaign's recorded `first_exact_basis` at $n=45$ and discharges the second stop condition of `bench/spec.json:next_research_campaign`; the VeriPB and $\mathrm{srg}(45,22,10,11)$ complementary-SAT tiers then remain the recorded residual fallback.

**Motivation and strategic context.**

Literature status (DS1 rev Apr 24 2026): $R(5,5)\in\{43,\dots,46\}$ — lower bound 43 (Exoo 1989), upper bound $46$ (Angeltveit–McKay, arXiv:2409.15709 v2). At $n=45$ the only published structural restriction is the degree window $\{20,\dots,24\}$. The separable $T$-family identity line is now CLOSED in-repo: $m=3$ ended in `M3_NO_CUT_IN_FROZEN_CONE`, $m=4$ in `M4_NO_CUT_IN_FROZEN_BASIS` — and the $m=4$ search proved every accepted certificate carries $\varepsilon=\zeta=0$, i.e. the $m=4$ row added zero cut strength over the $m=3$ basis. MR97 Lemma 2.3 exhausts all separable identities of degree $\le 6$, and $m\ge5$ is void for $(5,5)$-graphs. The only recorded unused identity in the campaign's first exact basis is Engström's mixed $K=4$ row, proved in arXiv:1002.4304 ("Theorem (Conjecture 1 in McKay and Radziszowski)"), extracted and machine-verified in `notes/engstrom_identity_extraction_2026-08-18.md`.

The identity, with $n=|V(G)|$ the global order, VERBATIM coefficients (extraction note §1, arxiv3.tex lines 79–104), is the universal per-vertex sum

$$
\sum_v \big[\,p_1(G_v^+)+p_2(G_v^-)+p_3(G_v^+,G_v^-)\,\big]=0,
$$

$$
\begin{aligned}
p_1(X)&=n(n-3)s(K_1,X)-(n^2+2n-6)s(K_1,X)^2+3n\,s(K_1,X)^3-2s(K_1,X)^4\\
&\quad+2(n^2+n-8)s(K_2,X)-12s(K_2,X)^2-12(n-1)s(K_1,X)s(K_2,X)\\
&\quad+12s(K_1,X)^2s(K_2,X)+72s(C_4,X)+12(n-2)s(K_3,X)+24s(K_{1,3},X)\\
&\quad+24s(P_4,X)+24s(T_{3,1},X)+12(n+2)s(P_3,X)-24s(K_1,X)s(P_3,X)+32s(T_{3,2},X);\\
p_2(Y)&=4s(K_2,Y)^2-12s(K_{1,3},Y)-8s(C_4,Y)-8s(T_{3,1},Y)-24s(T_{3,2},Y)+2(n-8)s(P_3,Y);\\
p_3(X,Y)&=4\,s(K_1,X)\,s(P_3,Y)-2(n-2)\,s(K_1,X)\,s(K_2,Y)+4\,s(K_1,X)^2\,s(K_2,Y),
\end{aligned}
$$

where $s(J,H)$ counts induced copies, $T_{3,1}$ is the paw, $T_{3,2}$ is the diamond, and every monomial has total argument-vertex count 4 (the $K=4$ case). Machine verification on record: residual zero on all 33,867 labeled graphs $n\le6$ plus 300 random graphs at $n\in\{7,\dots,10\}$. Unlike the $T$-family rows, this identity needs **no** $s(K_4,\cdot)$ counter at all and needs **no** $K_4$-free hypothesis — it holds for every graph, so the kernel has exactly one (universal) form.

Why this row is genuinely new machinery: its coordinate support ($C_4$, $K_{1,3}$, $P_3$, $P_4$ on both sides, plus three mixed product coordinates $d_v\cdot s(P_3,G_v^-)$, $d_v\cdot s(K_2,G_v^-)$, $d_v^2\cdot s(K_2,G_v^-)$) is disjoint from the frozen basis rows' support ($\{d,e,t,\mathrm{paw},\mathrm{diamond}\}$ on the $N$-side, $i_4$ on the $D$-side for $m=4$), and the $p_3$ cross terms make it non-separable, hence outside MR97 Lemma 2.3. Engström calls it "the first example of one outside [the easily described families]".

Measured expectations are deliberately muted (see the OPEN RISK section of this plan): the frozen cone's history says new rows may still add nothing. At $n=49$ the per-vertex Engström values on the four $(X,\bar Z)$ type combos are $\{0,144,288,432\}$ — all multiples of 144, non-negative, unique zero corner — so the row does NOT reproduce the MR97 Thm 3.2 contradiction there; it is a candidate CUT row at $n=45$, and the $n=49$ numbers serve as a kernel ground-truth replay fixture.

Fallback order recorded (executed only under separately approved designs): (1) VeriPB certificates for the existing gluing runs; (2) **this** tier — the Engström mixed identity, executed by the present plan; (3) complementary SAT on $\mathrm{srg}(45,22,10,11)$-type configurations. On a negative disposition the recorded residual order for the terminal record is VeriPB $\to$ $\mathrm{srg}(45,22,10,11)$ complementary SAT.

**Architecture:** First implement and exhaustively validate the universal mixed-identity kernel: all eleven raw induced 4-vertex motif counters (eight new; $i_4$, diamond, $K_3{+}K_1$, $K_4$ reused from the landed $m=4$ kernel), the corrected complement-pair map, the exact per-vertex row, and the residual, with the $n=49$ type-combo and the 656-graph $n=42$ ground-truth replays. Then stream every published-complete $R(4,5)$ catalog at orders 17–24 once, recording exact per-class histograms of fifteen aggregates, derive sound windows for missing classes (envelope LP for the eleven raw coordinates over the corrected R1–R5 relation system; degree-histogram aggregates for the 3-vertex-derived ones), and extend the frozen $n=45$ state cone — retaining the frozen $m=2/m=3$ rows $g$ and the $m=4$ row $h$ — with the mixed-row interval $F$. Search exact affine certificates over the enlarged basis with eight coefficients, and re-verify everything with a disjoint checker that recomputes the $F$ endpoints from its own kernel. Every certificate carries the $m=3$, $m=4$, and mixed rows.

**Tech Stack:** Python 3.14.3, standard-library `unittest`, integer bitsets, `fractions.Fraction`, SciPy 1.18.0 `linprog` for coefficient discovery only, validated graph6 catalogs, existing frozen modules `check_ramsey`, `subgraph_identities`, `m3_deficiency_cone`, `m4_subgraph_identities`, `m4_deficiency_cone`.

## Trust Roots

The six trust roots of the landed $m=4$ campaign carry over verbatim (recorded in `data/higher_identity_m4.json:trust_roots`; the first four are the $m=3$ roots):

1. Each selected order-17..23 file of `data/r45extreme` is the complete census of one fixed-edge $R(4,5)$ class (Angeltveit–McKay census, arXiv:2409.15709); not re-proved locally.
2. The 13 catalogs `r35_1.g6..r35_13.g6` are the complete censuses of $R(3,5)$ graphs (McKay data page); not re-proved locally.
3. `r45_24.g6` is the complete census of the 352,366 $R(4,5,24)$ graphs (arXiv:1703.08768); not re-proved locally.
4. `r55_42some.g6` holds 328 published $R(5,5,42)$ representatives; replay proves the identity on them, not that they exhaust $R(5,5,42)$.
5. The $m=4$ identity: locally re-derived, exhaustively verified $n\le6$, spot-verified $7$–$8$, replayed on all 656 $n{=}42$ graphs and complements and against the MR97 Thm 3.2 $n=49$ numbers (anchors MR97 Thm 2.2/3.2, https://users.cecs.anu.edu.au/~bdm/papers/r55.pdf).
6. The motif-window tables in `data/higher_identity_m4.json`: hash-pinned frozen catalogs plus the relation system R1–R5; consumed downstream only through the artifact.

New mixed-tier trust roots:

7. **The Engström $K=4$ mixed identity itself**: locally re-derived (module-docstring counting ledger, Task 1), exhaustively verified on all labeled graphs $n\le6$ (33,867 graphs), spot-verified $7$–$8$, replayed on all 656 $n{=}42$ graphs and complements, and replayed against the recorded $n=49$ type-combo values $\{0,144,288,432\}$ (unique zero corner) on the frozen order-24 census. Primary-source anchors: Engström, arXiv:1002.4304, "Theorem (Conjecture 1 in McKay and Radziszowski)" (arxiv3.tex lines 79–104); MR97 Conjecture 1, https://users.cecs.anu.edu.au/~bdm/papers/r55.pdf.
8. **The local motif-window tables** in `data/engstrom_identity.json`: derived only from hash-pinned frozen catalogs, the relation system R1–R5 (whose constants are unit-tested), and the frozen $m=3$ $q$-machinery; consumed downstream only through the artifact.

## Global Constraints

- Canonical design: `math/r55/bench/spec.json:next_research_campaign`; this plan executes the Engström mixed-tier slice of its `first_exact_basis` (fourth listed identity). The slot status is `M4_NO_CUT_IN_FROZEN_BASIS`; the basis-level stop condition stays OPEN until this tier is executed and recorded. This plan must not start before the $m=4$ campaign has landed (it owns the frozen artifacts and the module patterns this plan mirrors); it has.
- AI and SciPy remain outside the trust base; every accepted coefficient is rationalized and checked exactly.
- Reuse `check_ramsey` for graph6 parsing and `popcount`; `subgraph_identities` for `induced_subgraph`, `complement`, `edge_count`, `triangle_count`, `induced_p3_count`, `independent_triple_count`; `m4_subgraph_identities` for `diamond_count`, `independent_quad_count`, `induced_k3k1_count`, `k4_count`; `m4_deficiency_cone` for `affine_interval`, the R1–R5 relation rows, the relaxed-RHS enumeration, and the $g/h$ state endpoints. Do not create a second parser, a second triangle counter, or a second copy of any landed counter.
- The frozen artifacts `data/higher_identity_m3.json` and `data/higher_identity_m4.json` are read-only evidence; mixed-tier results go to a fresh artifact `data/engstrom_identity.json` at `schema_version: 3`.
- **Corrected-math rule (v3-typed).** The exact counter forms specified in Task 1 supersede every research-note shortcut. Six earlier candidate forms are provably WRONG and are pinned as detected mutations in the Task 1 tests: (i) $\mathrm{claw}=\sum_v\binom{d_v}{3}-\mathrm{paw}$ (drops $-2\cdot\mathrm{diamond}-4K_4$; the 4-vertex diamond gives $\sum\binom{d_v}{3}{=}\,2\neq0$); (ii) $K_{1,3}(\bar Z)=\mathrm{paw}(Z)$ (correct partner: $s(K_3{+}K_1,Z)$; witness $Z=K_3{+}K_1$); (iii) $C_4(\bar Z)=C_4(Z)$ ($C_4$ is NOT self-complementary; correct partner: $2K_2(Z)$; witness $Z=2K_2$); (iv) $s(P_3,\bar Z)=\sum_{uv\in E(Z)}(\nu-d_u-d_v)$ (missing $+c_{uv}$; witness $Z=K_3$ gives $-3\neq0$); (v) $K_2{+}2K_1=\sum_{uv\in E}\binom{\nu-|N[u]\cup N[v]|}{2}$ (overcounts by $2\cdot(2K_2)$; witness $Z=2K_2$); (vi) $\mathrm{paw}=\sum_v(d_v-2)t_v-4\cdot\mathrm{diamond}$ as a UNIVERSAL form (missing $-12K_4$; witness $Z=K_4$). On the swept $R(4,5)$ catalogs $K_4\equiv0$ by Ramseyhood — asserted during streaming, never assumed. Fast-versus-naive cross-validation on all labeled graphs $n\le6$ is the hard enforcement.
- v3 type pins: the cone `State` has exactly eleven integer fields `(d, a, b, deficiency, excess_balance, g_lo, g_hi, h_lo, h_hi, f_lo, f_hi)`; a `Certificate` carries `objective_id` plus exactly eight canonical reduced-fraction strings `(alpha, beta, gamma, delta, epsilon, zeta, eta, theta)`; the artifact is exactly `data/engstrom_identity.json` at schema 3 with the fifteen histogram keys frozen in Task 2.
- Missing inputs, count drift, identity mismatch, replay mismatch, or certificate mismatch must exit nonzero.
- Do not change the frozen benchmark engine, the candidate-3 artifacts, the frozen $m=3/m=4$ producers, or any frozen artifact; do not add dependencies.
- No $m=5$ $T$-family work (void), no Engström $K\ge5$ work (Conjectures 1.1/1.2 are open; out of scope), no SDP, no gluing engine, and no VeriPB or SAT implementation in this plan; the fallback order is recorded, not executed.
- No commit, push, or external mutation without separate authorization; this project tree is currently untracked by Git.
- A valid negative result is `MIXED_NO_CUT_IN_FROZEN_BASIS`, backed by exact primal witnesses for every instantiable route; never relabel a numerical trend, unresolved certification, or relaxation result as a theorem about all possible uses of the mixed identity. On that disposition the terminal record states additionally that the campaign's `first_exact_basis` ($m=2,3,4$ plus the mixed identity) is now exhausted, fulfilling the second stop condition of `next_research_campaign`, and records the residual fallback order (VeriPB $\to$ $\mathrm{srg}(45,22,10,11)$ complementary SAT) as the decision point for the next separately reviewed campaign.
- Do not claim $R(5,5)\le$ anything new.

---

### Task 1: Exact mixed-identity kernel

**Files:**
- Create: `math/r55/src/mixed_subgraph_identities.py`
- Create: `math/r55/tests/test_mixed_subgraph_identities.py`

**Interfaces:**
- Consumes: `check_ramsey.popcount`; `subgraph_identities.{complement, induced_subgraph, edge_count, triangle_count, induced_p3_count, independent_triple_count, motif3}`; `m4_subgraph_identities.{diamond_count, independent_quad_count, induced_k3k1_count, k4_count, motif4}`.
- Produces (all exact-integer; the eight new counters are the fast forms below, each cross-validated against naive 4-subset enumeration):

```python
MotifSweep = namedtuple(
    "MotifSweep",
    "order edges wedges triangles induced_p3 independent_triples "
    "e4 k2_2k1 two_k2 p3_k1 p4 claw k3k1 c4 paw diamond k4 pc3",
)

def c4_count(adj: list[int]) -> int            # induced cycles on 4 vertices
def claw_count(adj: list[int]) -> int          # induced K_{1,3}
def paw_count(adj: list[int]) -> int           # induced paw (T_{3,1})
def two_k2_count(adj) -> int                   # induced 2K2 (two disjoint edges)
def k2_2k1_count(adj) -> int                   # induced K2 + 2K1 (single edge)
def p3_k1_count(adj) -> int                    # induced P3 + K1
def p4_count(adj) -> int                       # induced paths on 4 vertices
def complement_p3_count(adj) -> int            # s(P3, complement(adj)); = pc3
def motif_sweep(adj: list[int]) -> MotifSweep  # one O(nu^2)-class pass, all fields
def mixed_vertex_row(n: int, d: int, x: MotifSweep, z: MotifSweep) -> int
def mixed_residual(adj: list[int]) -> int      # sum_v mixed_vertex_row(...)
```

`mixed_vertex_row(n, d, x, z)` computes $F(X,Y)$ for $X=G_v^+$ of order $d$ and $Z=\overline{G_v^-}$ of order $m=n-1-d$ (assert `x.order == d` and `x.order + z.order == n - 1`, mirroring the $m=4$ row's order-sum guard), applying the complement map internally so the Y-side polynomials are single-sourced:

```text
s(K2, Y) = C(m,2) - z.edges                      (exact)
s(P3, Y) = z.pc3                                  (window key "pc3")
s(C4, Y) = z.two_k2                               (window key "two_k2")
s(K1,3, Y) = z.k3k1                              (window key "k3k1")
s(paw, Y) = z.p3_k1                              (window key "p3_k1")
s(diamond, Y) = z.k2_2k1                        (window key "k2_2k1")
```

`mixed_residual(adj) -> int` sums `mixed_vertex_row(n, d_v, motif_sweep(N_v), motif_sweep(complement(D_v)))` over all $v$; there is exactly one form (the identity is universal; no $K_4$-free hypothesis is needed or used).

Fast closed forms (production path; `nu` = local order of the graph argument):

$$
\begin{aligned}
\mathrm{p}_3&= \text{(imported) }\ \textstyle\sum_v\binom{d_v}{2}-3t,\\
\mathrm{paw}&=\textstyle\sum_v (d_v-2)\,t_v-4\,\mathrm{diamond}-12\,K_4,\qquad t_v=\#\text{triangles through }v,\\
\mathrm{c}_4&=\tfrac12\textstyle\sum_{uv\notin E}\big[\tbinom{c_{uv}}{2}-e(G[\mathrm{cn}(u,v)])\big],\qquad c_{uv}=|\mathrm{cn}(u,v)|,\ \ \mathrm{cn}=\mathrm{common\ neighbors},\\

\emph{c4 correction.} The earlier $\ell_2=\tfrac12\#\{uv\notin E:c_{uv}=2\}$ threshold form is WRONG and is superseded: a non-edge with $c_{uv}>2$ common neighbors participates in $\binom{c_{uv}}{2}-e(G[\mathrm{cn}])$ distinct C4s (e.g.\ $K_{2,3}$ has $c_4=3$ but $\ell_2=3$, giving $3\neq(3-0)/2=1$). Each non-adjacent common-neighbor pair completes one C4 with the non-edge; adjacent ones close diamonds (which have only one non-edge, contributing nothing). The identity $\sum_{uv\notin E}[\binom{c_{uv}}{2}-e(G[\mathrm{cn}])]=2\,c_4$ is the correct per-vertex incidence (enforced exhaustively $n\le6$ in the landed tests).\\
\mathrm{claw}&=\textstyle\sum_v\binom{d_v}{3}-\mathrm{paw}-2\,\mathrm{diamond}-4\,K_4,\\
2K_2&=\tfrac12\textstyle\sum_{uv\in E} e\big(G[\mathrm{nonN}(u)\cap\mathrm{nonN}(v)]\big),\\
K_2{+}2K_1&=\textstyle\sum_{uv\in E}\binom{\nu-|N[u]\cup N[v]|}{2}-2\cdot(2K_2),\\
P_3{+}K_1&=\textstyle\sum_x \mathrm{induced\_p3\_count}\big(G-N[x]\big),\qquad \text{(}O(\nu^3)\text{, bitset subgraph per }x\text{)},\\
P_4&=\tbinom{\nu}{4}-\text{(sum of the other ten induced 4-vertex classes)},\\
\mathrm{pc}_3&=\text{(form B, sweep path)}\ =\nu\,e(\text{input})... \ \text{see the three equivalent forms below}.
\end{aligned}
$$

`pc3` is pinned by THREE provably equivalent closed forms, all enforced equal exhaustively: (A) $\sum_{uv\in E(Z)}(\nu-d_u-d_v+c_{uv}(Z))$; (B) $e(Z)(\nu-2)-2\,s(P_3,Z)-3\,t(Z)$; (C) $\sum_v\binom{\nu-1-d_v}{2}-3\,i_3(Z)$ with $i_3=\binom{\nu}{3}-e(\nu-2)+\sum_v\binom{d_v}{2}-t$. Form (B) is the cheapest and is the sweep path; (A) is Main's section-5 line. Dropping the $+c_{uv}$ term (candidate form iv) must be detected on witness $Z=K_3$ (residual $-3\neq0$).

**Universal complement-pair map (extraction note §5, corrected in-session 2026-08-18/19; exhaustive n<=6 tests pin each line):**

$$
\begin{aligned}
s(P_3,\bar Z)&=\mathrm{pc}_3(Z), & C_4(\bar Z)&=2K_2(Z), & K_{1,3}(\bar Z)&=s(K_3{+}K_1,Z),\\
\mathrm{paw}(\bar Z)&=s(P_3{+}K_1,Z), & \mathrm{diamond}(\bar Z)&=s(K_2{+}2K_1,Z), & K_4(\bar Z)&=s(E_4,Z)=i_4(Z),\\
P_4(\bar Z)&=P_4(Z), & e(\bar Z)&=\binom{\nu}{2}-e(Z).
\end{aligned}
$$

No complement graph is ever materialized for window purposes; the kernel materializes complements only in tests.

**Module docstring counting ledger (required).** The docstring must contain the independent derivation in the ledger form of these verified anchor lemmas: (L1) $\sum_v s(J,N_v)=\sum_{S\cong J}|\mathrm{cn}(S)|$ and $\sum_v s(J,D_v)=\sum_{S\cong J}|\mathrm{cnn}(S)|$; (L2) cross products $\sum_v d_v^{k}s(J,D_v)=\sum_{S\cong J}\sum_{v\in\mathrm{cnn}(S)}d_v^{k}$, $k\in\{1,2\}$ (and dually on the $N$-side); (L3) same-side products expand as $\sum_{S_1\cong J_1}\sum_{S_2\cong J_2}|\cdots|$ regrouped by the induced type on $S_1\cup S_2$, which stays within four vertices by the $K=4$ bound; (L4) the 25 monomials (16 in $p_1$, 6 in $p_2$, 3 in $p_3$) are regrouped so every degree-1 normal-form coefficient of Engström's Theorem 2.18/Corollary 2.19 vanishes, $\sum_v F_v=0$ closes pattern-by-pattern; (L5) the exhaustive labeled-graph test is the second, mechanically independent check. Also record the three pc3 forms and the pinned witness graphs of the corrected-math rule.

- [ ] **Step 1: Write the exhaustive failing tests**

Create `test_mixed_subgraph_identities.py` with the shared `graph_from_mask` helper, a naive reference 4-subset motif counter (test-only), and these observable contracts:

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
from m4_subgraph_identities import (
    diamond_count, independent_quad_count, induced_k3k1_count, k4_count,
)
from mixed_subgraph_identities import (  # expected to fail before Task 1 impl
    MotifSweep, c4_count, claw_count, complement_p3_count, k2_2k1_count,
    mixed_residual, mixed_vertex_row, motif_sweep, p3_k1_count, p4_count,
    paw_count, two_k2_count,
)


def graph_from_mask(n, mask):
    adj = [0] * n
    for bit, (u, v) in enumerate(itertools.combinations(range(n), 2)):
        if (mask >> bit) & 1:
            adj[u] |= 1 << v
            adj[v] |= 1 << u
    return adj


class MixedIdentityTests(unittest.TestCase):
    def test_fast_counters_match_naive_through_six_vertices(self):
        for n in range(2, 7):
            for mask in range(1 << (n * (n - 1) // 2)):
                adj = graph_from_mask(n, mask)
                self.assertEqual(c4_count(adj), naive_c4(adj))
                self.assertEqual(claw_count(adj), naive_claw(adj))
                self.assertEqual(paw_count(adj), naive_paw(adj))
                self.assertEqual(two_k2_count(adj), naive_two_k2(adj))
                self.assertEqual(k2_2k1_count(adj), naive_k2_2k1(adj))
                self.assertEqual(p3_k1_count(adj), naive_p3_k1(adj))
                self.assertEqual(p4_count(adj), naive_p4(adj))
                self.assertEqual(complement_p3_count(adj), naive_p3(complement(adj)))

    def test_complement_pair_map_exhaustive_through_six_vertices(self):
        for n in range(2, 7):
            for mask in range(1 << (n * (n - 1) // 2)):
                adj = graph_from_mask(n, mask)
                cj = complement(adj)
                self.assertEqual(naive_p3(cj),  complement_p3_count(adj))
                self.assertEqual(naive_c4(cj),  two_k2_count(adj))
                self.assertEqual(naive_claw(cj), induced_k3k1_count(adj))
                self.assertEqual(naive_paw(cj),  p3_k1_count(adj))
                self.assertEqual(naive_diamond(cj), k2_2k1_count(adj))
                self.assertEqual(naive_k4(cj),  independent_quad_count(adj))
                self.assertEqual(naive_p4(cj),   p4_count(adj))

    def test_pc3_three_closed_forms_agree_exhaustive(self):
        for n in range(2, 7):
            for mask in range(1 << (n * (n - 1) // 2)):
                adj = graph_from_mask(n, mask)
                a = pc3_form_a(adj)   # sum over edges of (nu - du - dv + cuv)
                b = edge_count(adj) * (len(adj) - 2) - 2 * induced_p3_count(adj) - 3 * triangle_count(adj)
                c = pc3_form_c(adj)
                self.assertEqual(a, b); self.assertEqual(b, c)
        # dropped-c_uv mutation must NOT agree: pinned witness K3
        k3 = graph_from_mask(3, 0b111)
        self.assertEqual(pc3_form_a(k3), 0)
        self.assertEqual(pc3_form_a_without_cuv(k3), -3)

    def test_mixed_residual_exhaustive_through_six_vertices(self):
        total = 0
        for n in range(1, 7):
            for mask in range(1 << (n * (n - 1) // 2)):
                self.assertEqual(mixed_residual(graph_from_mask(n, mask)), 0)
                total += 1
        self.assertEqual(total, 33867)

    def test_spot_random_seven_eight_vertices(self):
        rng = random.Random(20260818)
        for _ in range(40):
            n = rng.choice([7, 8])
            adj = [0] * n
            for u, v in itertools.combinations(range(n), 2):
                if rng.random() < rng.choice([0.3, 0.5, 0.7]):
                    adj[u] |= 1 << v
                    adj[v] |= 1 << u
            self.assertEqual(mixed_residual(adj), 0)

    def test_mutated_coefficients_are_detected(self):
        for mutation in ("c4_coef", "cross_coef", "p2_c4_coef"):
            worst = max(
                abs(mutated_residual(graph_from_mask(n, mask), mutation))
                for n in range(4, 7)
                for mask in range(1 << (n * (n - 1) // 2))
            )
            self.assertGreater(worst, 0, mutation)

    def test_closed_form_corrections_are_detected(self):
        # pinned witnesses of the corrected-math rule; each wrong form yields
        # a nonzero discrepancy on an exact small graph
        diamond = graph_from_mask(4, 0b111110)      # K4 minus edge (0,1)
        self.assertEqual(sum(_c3(deg) for deg in (3, 3, 2, 2)), 2)
        self.assertEqual(claw_count(diamond), 0)     # = 2 - paw(0) - 2*dia(1) - 4*K4(0)
        two_k2 = graph_from_mask(4, 0b010010)       # edges: (0,2), (1,3)
        self.assertEqual(k2_2k1_count(two_k2), 0)   # bare C(2,2)*2 = 2, minus 2*(2K2)=2
        k4 = graph_from_mask(4, (1 << 6) - 1)
        self.assertEqual(paw_count(k4), 0)           # = 12 - 0*dia - 12*K4

The false invariant "sum_x e(G - N[x]) = 4*(2K2) + 2*(K2+2K1)" and the
"l2 == 2*C4 + diamond" audit have been REMOVED from the landed test list
(falsified: for $P_3+K_1$, sum_x e(G-N[x]) = 2 with a 4-element
closed-neighborhood complement of size 3, not 0); the landed tests assert
only the correct per-vertex incidence sum_{uv not in E}[C(c_uv,2) -
e(G[cn])] = 2 c4.

    def test_all_656_known_r55_graphs_have_zero_mixed_residual(self):
        published = []
        with (ROOT / "data" / "r55_42some.g6").open(encoding="ascii") as fh:
            for line in fh:
                if line.strip():
                    _, adj = parse_graph6_line(line)
                    published.append(adj)
        self.assertEqual(len(published), 328)
        tested = 0
        for adj in published:
            self.assertEqual(mixed_residual(adj), 0)
            self.assertEqual(mixed_residual(complement(adj)), 0)
            tested += 2
        self.assertEqual(tested, 656)

    def test_n49_type_combo_replay(self):
        graphs = []
        with (ROOT / "data" / "r45_24.g6").open(encoding="ascii") as fh:
            for line in fh:
                if not line.strip():
                    continue
                _, adj = parse_graph6_line(line)
                if sum(row.bit_count() for row in adj) // 2 == 132:
                    graphs.append(adj)
        self.assertEqual(len(graphs), 2)
        for x in graphs:                       # 11-regular; Thm 3.2 anchor
            sw = motif_sweep(x)
            self.assertEqual(sw.triangles, 176)
            self.assertEqual(sw.diamond, 792)
            self.assertEqual(sw.paw, 1584)     # new deliberate pin on X1/X2
            self.assertEqual(sw.claw, 792)      # = sum C(11,3)*24 - paw - 2*diamond
            self.assertEqual(sw.k3k1, 528)
        values = sorted(
            mixed_vertex_row(49, 24, motif_sweep(xi), motif_sweep(zj))
            for xi in graphs for zj in graphs
        )
        self.assertEqual(values, [0, 144, 288, 432])
        self.assertEqual(sum(v == 0 for v in values), 1)   # unique zero corner
```

`naive_*` (4-subset enumeration), `pc3_form_a` (with and without the `c_uv` term), `pc3_form_b`, `pc3_form_c`, `_c3`, and `mutated_residual` (mutation selectors: `"c4_coef"` perturbs the $72s(C_4,X)$ coefficient to 71; `"cross_coef"` perturbs $-2(n-2)$ to $-2(n-1)$ in $p_3$; `"p2_c4_coef"` perturbs $-8s(C_4,Y)$ to $-7$) live in the test file.

- [ ] **Step 2: Run the test and observe the intended failure**

```bash
cd /Users/jinleic/jinleic-workspace/math
./.venv/bin/python -m unittest r55/tests/test_mixed_subgraph_identities.py -v
```

Expected: import failure for `mixed_subgraph_identities`.

- [ ] **Step 3: Implement the exact kernel**

Implementation shape in `mixed_subgraph_identities.py` — use the landed fast forms from `m4_subgraph_identities` untouched; the eight new counters are:

The landed kernel supersedes the draft pseudocode on four points (the draft's
own exhaustive tests caught each; the landed forms are given in
`mixed_subgraph_identities.py`):

1. `paw`'s `t_v = _pairs_inside(adj, adj[v])` (edges inside N_v), not a
   two-pointer triangle-enumeration pass.  The draft two-pointer loop
   contributed each triangle through v to `tri_through[v]` twice and
   omitted any pass that counted vertex-degree-2 triangles.
2. `c4` uses the corrected $\frac12\sum_{uv\notin E}[\binom{c_{uv}}2-
   e(G[\mathrm{cn}])]$ form above, not the $\ell_2$ threshold.
3. `k2_2k1`'s closed neighborhood $N[u]\cup N[v]$ is `popcount(adj[u] |
   adj[v])` with NO `+2`: $u$ and $v$ are both already represented in the
   union bitmasks (the draft's `+2` double-counted them; its own
   `k2_2k1(2K2)=0` pin caught it: bare sum $2$, and $2-2\cdot 1=0$).
4. `complement_p3_count(adj) = m3.edges*(n-2) - 2*m3.induced_p3 -
   3*m3.triangles` computed from a single `motif3` pass; calling
   `motif_sweep` here would be circular.

```python
def paw_count(adj):
    n = len(adj)
    tri_through = [_pairs_inside(adj, adj[v]) for v in range(n)]
    return (sum((popcount(adj[v]) - 2) * tv for v, tv in enumerate(tri_through))
            - 4 * diamond_count(adj) - 12 * k4_count(adj))

def claw_count(adj):
    c3 = sum(popcount(r) * (popcount(r) - 1) * (popcount(r) - 2) // 6 for r in adj)
    return c3 - paw_count(adj) - 2 * diamond_count(adj) - 4 * k4_count(adj)

def c4_count(adj):
    n = len(adj); total = 0
    for u in range(n):
        for v in range(u + 1, n):
            if (adj[u] >> v) & 1: continue
            cn = adj[u] & adj[v]; c = popcount(cn)
            total += c * (c - 1) // 2 - _pairs_inside(adj, cn)
    return total // 2

def two_k2_count(adj):
    n = len(adj); t = 0; full = (1 << n) - 1
    for u in range(n):
        for v in range(u + 1, n):
            if not ((adj[u] >> v) & 1): continue
            mask = full & ~adj[u] & ~adj[v] & ~((1 << u) | (1 << v))
            t += _pairs_inside(adj, mask)
    return t // 2

def k2_2k1_count(adj):
    n = len(adj); s = 0
    for u in range(n):
        for v in range(u + 1, n):
            if not ((adj[u] >> v) & 1): continue
            closed = popcount(adj[u] | adj[v])   # |N[u] u N[v]|, no +2
            outside = n - closed
            s += outside * (outside - 1) // 2
    return s - 2 * two_k2_count(adj)

def p3_k1_count(adj):
    n = len(adj); t = 0
    for x in range(n):
        keep = [u for u in range(n) if u != x and not ((adj[x] >> u) & 1)]
        t += induced_p3_count(induced_subgraph(adj, keep))
    return t

def p4_count(adj):
    n = len(adj)
    return _c4(n) - (independent_quad_count(adj) + k2_2k1_count(adj)
                      + two_k2_count(adj) + p3_k1_count(adj)
                      + induced_k3k1_count(adj) + claw_count(adj)
                      + c4_count(adj) + paw_count(adj)
                      + diamond_count(adj) + k4_count(adj))

def complement_p3_count(adj):
    n = len(adj); m3 = motif3(adj)
    return m3.edges * (n - 2) - 2 * m3.induced_p3 - 3 * m3.triangles

def motif_sweep(adj):
    m3 = motif3(adj)
    m4 = motif4(adj)
    return MotifSweep(
        order=len(adj), edges=m3.edges,
        wedges=sum(popcount(r) * (popcount(r) - 1) // 2 for r in adj),
        triangles=m3.triangles, induced_p3=m3.induced_p3,
        independent_triples=m3.independent_triples,
        e4=m4.independent_quads, k2_2k1=k2_2k1_count(adj),
        two_k2=two_k2_count(adj), p3_k1=p3_k1_count(adj), p4=p4_count(adj),
        claw=claw_count(adj), k3k1=m4.k3k1, c4=c4_count(adj),
        paw=paw_count(adj), diamond=m4.diamonds, k4=m4.cliques,
        pc3=complement_p3_count(adj),
    )
```

The module docstring must state one counting derivation per counter. `mixed_vertex_row` is the single source of the coefficient layout:

```python
def mixed_vertex_row(n, d, x, z):
    assert x.order == d and x.order + z.order == n - 1
    ex = x.edges                    # s(K2, X) = e(G_v^+)
    ey = _c2(n - 1 - d) - z.edges  # s(K2, Y) = e(G_v^-), Y = complement(Z)
    exact = (n * (n - 3) * d - (n * n + 2 * n - 6) * d * d + 3 * n * d * d * d
             - 2 * d ** 4 + 2 * (n * n + n - 8) * ex - 12 * ex * ex
             - 12 * (n - 1) * d * ex + 12 * d * d * ex + 4 * ey * ey
             - 2 * (n - 2) * d * ey + 4 * d * d * ey)
    return exact + (72 * x.c4 + 12 * (n - 2) * x.triangles + 24 * x.claw
                   + 24 * x.p4 + 24 * x.paw
                   + (12 * (n + 2) - 24 * d) * x.induced_p3 + 32 * x.diamond)
                   + (-12 * z.k3k1 - 8 * z.two_k2 - 8 * z.p3_k1
                      - 24 * z.k2_2k1 + (2 * (n - 8) + 4 * d) * z.pc3)

def mixed_residual(adj):
    n = len(adj)
    return sum(mixed_vertex_row(n, d_v, motif_sweep(N_v),
                               motif_sweep(complement(D_v)))
               for v in range(n))
```

A `ValueError` guards `x.order + z.order != n - 1`.

- [ ] **Step 4: Run exhaustive tests and the existing graph-parser selftests**

```bash
cd /Users/jinleic/jinleic-workspace/math
./.venv/bin/python -m unittest r55/tests/test_mixed_subgraph_identities.py -v
./.venv/bin/python r55/src/check_ramsey.py --selftest
./.venv/bin/python -m unittest r55/tests/test_subgraph_identities.py -v
./.venv/bin/python -m unittest r55/tests/test_m4_subgraph_identities.py -v
```

Expected: all new tests `OK` (including the 656-graph replay, the $n=49$ type-combo replay with the unique zero corner, and every pinned corrected-math witness); existing suites unchanged and `selftest OK`.

---

### Task 2: Exact motif windows and the extended $n=45$ state cone

**Files:**
- Create: `math/r55/src/mixed_deficiency_cone.py`
- Create: `math/r55/tests/test_mixed_deficiency_cone.py`
- Create at execution time: `math/r55/data/engstrom_identity.json` (schema 3; fresh file; the frozen `higher_identity_m3.json` and `higher_identity_m4.json` stay read-only)

**Interfaces:**
- Consumes: Task 1 `motif_sweep`; `m4_deficiency_cone.{affine_interval, relation_rows, relaxed_rhs_intervals, enumerate_basic_solutions, catalog_motif_windows, outer_motif_window, m4_state_interval, _solve_square_int}` (frozen and landed — one owner per fact); `m3_deficiency_cone.{load_edge_bounds, r35_edge_windows, degree_histograms, q_outer_interval, budget2, degree_bounds}`; `data/structural_tables.json`, `data/VALIDATION.json`, `data/VALIDATION_extreme.json`, the 36 validated edge-extremal `r45extreme/*.g6` catalogs at orders 17–23, the complete `r45_24.g6` census, and the frozen artifacts `higher_identity_m3.json`, `higher_identity_m4.json` (read-only; m4 states used for the stability assert).
- Produces:

```python
def stream_catalog_motifs(path) -> dict[tuple[int, int], dict[str, dict[int, int]]]
def catalog_motif_windows(records) -> dict[tuple[int, int], dict[str, tuple[int, int]]]
def outer_motif_window(order, edges, key) -> tuple[int, int]   # all 15 keys
def mixed_state_interval(d, a, b, motif_windows=None, motif_sources=None) -> State
def build_states(motif_windows=None, motif_sources=None) -> list[State]
def state_record(state, motif_windows=None, motif_sources=None) -> dict
```

The fifteen swept aggregates per graph class record exactly these JSON keys (ASCII names, frozen; 11 raw 4-motifs + 3 derived 3-motif aggregates + pc3):

```text
q, t, p3, pc3, i4, k4, diamond, k3k1, c4, claw, paw, p4, two_k2, k2_2k1, p3_k1
```

with `q = wedges - t` (the frozen m3 aggregate), `p3 = induced P3`, `pc3 = s(P3, complement)`, `i4` = E4, `k4` = K4, `diamond` = K4-e, `k3k1` = K3+K1, `c4` = C4, `claw` = K1,3, `paw` = T31, `p4` = P4, `two_k2` = 2K2, `k2_2k1` = K2+2K1, `p3_k1` = P3+K1.

Per-graph budget (the whole-task cost driver; `nu` <= 24, all O(nu^2)-class bit-op passes):

| Piece | Operations |
|---|---|
| degrees, wedges, e, t, i3, p3 | existing motif3/motif4 passes (reused) |
| diamond, k3k1, i4, k4 | existing m4 fast counters (reused) |
| c4 (corrected form) | per non-edge: AND for the common-neighbor mask, popcount, plus the edge count inside that mask (one `_pairs_inside` per non-edge) => Nu^2*(1 + c) row-AND/popcount pairs; same cost class, budget unchanged |
| paw, claw | per vertex triangle pass (e AND+popcount) + bookkeeping |
| two_k2 | per edge: 2 ANDs + popcount + row-ANDs => E*(2 + cbar) |
| k2_2k1 | per edge: 2 ANDs on closed-neighborhood union + popcount => E*3 |
| p3_k1 | per vertex: induced subgraph + p3 pass => nu^3 bit-ops (~1.4e4 at nu=24) |
| p4, pc3 | bookkeeping (C(nu,4) minus others; closed form B) |

Same cost class as the landed m4 kernel. Measured anchors (Apple M3 Ultra, m4 campaign): 0.287 ms/graph (order 19, e=90), 0.421 ms/graph (order 23, e=119), 0.507 ms/graph (order 24); the eight new counters add a constant factor well under 2x in the same class. Catalog total 8,500,211 graphs = 8,147,845 in the 36 validated order-17..23 files + 352,366 in `r45_24.g6`; the m4 sweep measured **2310 s** single-process. Fail-closed budget: **7200 s single-process**. Streaming discipline mirrors m4: one graph6 line at a time, fold into per-(order, edges) histograms, retain no graph, hash every file first, never hold more than one adjacency list, and assert `k4 == 0` on every swept graph (Ramseyhood; assert, never assume). Multiprocessing is explicitly unnecessary; do not add it.

**Missing-class envelope policy.** For a missing class the sound windows for the eleven raw 4-motifs come from the certified envelope LP over the SAME frozen R1–R5 relation system (constants unit-tested), exactly as landed in the m4 campaign (same corrected constants: R3 in the partial-derivative form collapsing to W; R4 is the degree-histogram form; R5 carries triangle multiplicities tau=(1,1,2,4)). Reuse `m4_deficiency_cone.{relation_rows, relaxed_rhs_intervals, _solve_square_int}` verbatim — do not copy. The one generalization: extrema are requested for **all eleven** raw coordinates (not the m4 pick set (9,6,0)); enumerate the feasible basic solutions once per class in `Fraction` arithmetic (with R1 and R2 exact; the landed system is 11 nonnegativity + 2 equality (R1, R2) + 3 two-sided (R3, R4, R5, i.e. 6 one-sided rows), 17 inequality rows total, vertices arising by activating 9 of the 17) and certify min/max per coordinate. SciPy output is a pruning hint only; every accepted endpoint is the exact min/max over all verified basic solutions, else the class falls back to the trivially sound `[0, C(nu,4)]` and is recorded in `lp_fallback_classes()`. The 3-vertex-derived aggregates use degree-histogram aggregates from the frozen m3 enumeration: `q` and `w = sum C(d,2)` outer intervals per class; `t = w - q`; `p3 = q - 2t`; `pc3 = e(nu - 2) - 2*p3 - 3*t` — all sign-aware affine compositions over the SAME feasible histogram family already enumerated by the imported `q_outer_interval` machinery (never wider than the true joint support; this matches the m4 t-window policy exactly). Present classes always use the exact per-graph swept histograms (never recombined marginals — `pc3`, `p3`, `c4`, `claw`, `paw` are nonlinear in the joint within a class).

**The new F interval per state.** Each state $(d,a,b)$ reuses the frozen $g$ and $h$ endpoints from `m4_state_interval` (import, do not copy) and adds the mixed-row interval $F_v\in[F_{lo},F_{hi}]$ over the two strata $x$ ($R(4,5,d)$, $e_x=E_{45}(d)-a$) and $z=\overline{D_v}$ ($R(4,5,m)$, $e_z=E_{45}(m)-b$, $m=44-d$). At $n=45$ (v3-typed constants):

$$
F_v = \underbrace{1890d-2109d^2+135d^3-2d^4+4124e_x-12e_x^2-528\,d\,e_x+12d^2e_x+4e_y^2-86\,d\,e_y+4d^2e_y}_{\text{state-exact}}
+516\,t(x)+72\,c_4(x)+24\,\mathrm{claw}(x)+24\,p_4(x)+24\,\mathrm{paw}(x)+(564-24d)\,p_3(x)+32\,\mathrm{diamond}(x)
+(74+4d)\,\mathrm{pc}_3(z)-12\,k_3k_1(z)-8\cdot2K_2(z)-8\,(P_3{+}K_1)(z)-24\,(K_2{+}2K_1)(z),
$$

with $e_y=\tbinom{44-d}{2}-e_z$. Every windowed term enters through `affine_interval(coef, window)` (sign-aware; the $p_3$ coefficient $564-24d\in\{84,60,36,12,-12\}$ flips sign at $d=24$, exactly as the m4 $t$ coefficient did). `f_lo <= f_hi` for every state, and $\sum_v F_v=0$ gives the new global pair $\sum_v F_{lo}\le0\le\sum_v F_{hi}$.

**Cone-stability invariant:** `build_states` re-enumerates the $(d,a,b)$ domain with the frozen m3 feasibility pruning and MUST reproduce the frozen m4 artifact's 3,215 state records on the projection `(d, a, b, g_lo, g_hi, h_lo, h_hi)` plus interval sources, keyed by `(d, a, b)` — a hard assertion, never a silent re-import. Missing $R(4,5)$ catalogs are never treated as infeasibility.

- [ ] **Step 1: Write failing envelope tests**

Tests must cover (naive reference counters allowed as test helpers):

1. `motif_sweep` equals the naive 18-field vector on 250 consecutive graphs of `r45extreme/r4517.77.g6` and on the two 132-edge order-24 graphs;
2. `affine_interval(-3, (10, 20)) == (-60, -30)` (imported helper, re-pinned);
3. every relation R1–R5 holds with the exact constants on 60 seeded random graphs ($n\le7$), and the imported `relation_rows`/`relaxed_rhs_intervals` agree with the landed m4 module outputs;
4. exact catalog windows override outer windows, for each of the fifteen keys;
5. for the three fixture stratum pairs of the m4 tests — $(20,100)$x$(24,132)$, $(21,107)$x$(23,122)$, $(22,114)$x$(22,114)$ — the exact $F$ value of every catalog-graph pair (first 20 of each side) lies inside `mixed_state_interval(...).f_lo..f_hi` (exact F computed via Task 1 `mixed_vertex_row`), and the $g,h$ intervals still contain the m3/m4 expressions;
6. `outer_motif_window` returns windows containing every catalog graph's value for each of the fifteen keys on at least three present classes (audit direction);
7. the state count matches the m4 enumeration domain ($d\in\{20,\dots,24\}$, full $(a,b)$ boxes; 3,215 states) AND every state satisfies `f_lo <= f_hi`, `g_lo <= g_hi`, `h_lo <= h_hi`;
8. the 656-graph replay and the $n=49$ type-combo replay both re-run through the Task 1 kernel.

- [ ] **Step 2: Run tests and observe missing-interface failures**

```bash
cd /Users/jinleic/jinleic-workspace/math
./.venv/bin/python -m unittest r55/tests/test_mixed_deficiency_cone.py -v
```

Expected: import failures for the cone interfaces.

- [ ] **Step 3: Implement streaming windows and the extended State**

```python
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
```

`build_states` imports the m4 machinery for the $g,h$ endpoints (never recomputes them from a second copy) and adds the $F$ endpoints from the twelve windowed terms of the two strata. It may omit a pair only when outer degree-histogram constraints prove one side infeasible, exactly as in the m4 plan; missing $R(4,5)$ catalogs are never treated as infeasibility. The v3 artifact writes `searches: []`, `disposition: "MIXED_ANALYZED"` until Task 3.

- [ ] **Step 4: Add fail-closed input provenance and the v3 artifact**

The CLI must:

1. hash every input file with SHA-256 (including the frozen m3 and m4 artifacts);
2. load edge windows from the validated complete `r35_k.g6` records (reuse `r35_edge_windows`);
3. assert all 328 lines plus complements give 656 zero mixed residuals;
4. replay the $n=49$ type-combo ground truth (two 132-edge graphs; four combos; multiset `[0, 144, 288, 432]`; exactly one zero corner; pins on paw=1584, claw=792);
5. stream every validated order-17..23 file and the complete order-24 census, one line at a time, asserting `k4 == 0` per graph and recording exact per-class histograms of the fifteen aggregates;
6. verify every observed value of every aggregate lies inside its `outer_motif_window` (envelope self-audit for all fifteen keys);
7. use exact catalog windows for present classes, envelope/LP or degree-histogram windows for missing classes, trivial windows only on certification failure;
8. assert the rebuilt state projection `(d, a, b, g_lo, g_hi, h_lo, h_hi)` equals the frozen m4 artifact's 3,215 records;
9. write JSON only after all checks pass, sorted keys, deterministic lists, atomic replace via temporary sibling.

The v3 JSON schema is fixed before search:

```text
schema_version: 3
campaign_id: higher_order_identity_positive_deficiency_mixed
inputs: sorted records {relative_path, sha256, bytes, graph_count}
trust_roots: the eight statements of this plan's Trust Roots section
replay: {published_graphs: 328, complements: 328, residual_zero: 656}
n49_replay: {extremal_graphs: 2, combos: 4, row_values: [0, 144, 288, 432],
             zero_corners: 1}
r35_edge_windows: sorted records {order, edge_min, edge_max}
catalog_motif_histograms: sorted records
  {order, edges, source, graph_count,
   histograms: {q: [...], t: [...], p3: [...], pc3: [...], i4: [...],
                k4: [...], diamond: [...], k3k1: [...], c4: [...],
                claw: [...], paw: [...], p4: [...], two_k2: [...],
                k2_2k1: [...], p3_k1: [...]}}
catalog_motif_windows: sorted records {order, edges, source, windows}
outer_motif_windows: sorted records {order, edges, method, windows}
states: sorted records
  {d, a, b, deficiency, excess_balance, g_lo, g_hi, h_lo, h_hi, f_lo, f_hi,
   x_interval_source, y_interval_source, x_motif_source, y_motif_source}
searches: [] until Task 3
disposition: MIXED_ANALYZED until Task 3
```

No wall-clock timestamp or absolute path enters the canonical JSON.

- [ ] **Step 5: Run focused tests and full analysis**

```bash
cd /Users/jinleic/jinleic-workspace/math
./.venv/bin/python -m unittest r55/tests/test_mixed_deficiency_cone.py -v
./.venv/bin/python r55/src/mixed_deficiency_cone.py \
  --output r55/data/engstrom_identity.json
```

Expected terminal contract:

```text
MIXED REPLAY: 656/656 residual zero
MIXED N49 REPLAY: extremal=2 combos=4 values=[0,144,288,432] zero_corners=1
MIXED STREAM: 8500211 graphs in <seconds> (budget 7200)
MIXED ENVELOPES: catalog violations=0 outer-window audits=ALL
MIXED STATES: 3215 states (asserted equal to frozen m4 projection)
WROTE r55/data/engstrom_identity.json
```

---

### Task 3: Exact affine certificates over the balance/g/h/F basis

**Files:**
- Create: `math/r55/src/search_mixed_cuts.py`
- Create: `math/r55/tests/test_search_mixed_cuts.py`
- Modify at execution time: `math/r55/data/engstrom_identity.json`

**Interfaces:**
- Consumes: the Task 2 v3 artifact (`States` = eleven-field frozen dataclass `(d, a, b, deficiency, excess_balance, g_lo, g_hi, h_lo, h_hi, f_lo, f_hi)`) and the frozen objective registry (identical routes and acceptance edges to m3/m4; NO new routes are registered in this tier — the fifteen Task 2 motif keys enter the search exclusively through the already-built per-state `f_lo/f_hi` interval endpoints, never as search features).
- Produces an immutable `Certificate` with fields `objective_id`, `alpha`, `beta`, `gamma`, `delta`, `epsilon`, `zeta`, `eta`, `theta`; every coefficient is a canonical reduced-fraction string.
- Produces an immutable `PrimalWitness` with `objective_id` and sorted `(state_index, weight)` rational pairs.
- Produces `exact_affine_bound(states, objective_id) -> Certificate`, `verify_certificate(states, certificate) -> Fraction`, `verify_primal_witness(states, witness) -> Fraction`, `verify_search_record(states, record) -> str`.

The eight-coefficient algebra. Every certificate satisfies gamma, delta, epsilon, zeta, eta, theta >= 0 (alpha, beta free) and, for every state,

    objective(state) <= alpha + beta*balance(state)
                      + gamma*g_lo(state) - delta*g_hi(state)
                      + epsilon*h_lo(state) - zeta*h_hi(state)
                      + eta*f_lo(state) - theta*f_hi(state).

The sign gates on the new pair are eta >= 0 multiplying f_lo and theta >= 0 multiplying f_hi — the same lo-positive/hi-negative pairing as (gamma, delta) and (epsilon, zeta). Aggregating over the 45 vertices uses sum balance = 0, sum g_lo <= 0 <= sum g_hi, sum h_lo <= 0 <= sum h_hi, and the mixed row sum f_lo <= 0 <= sum f_hi, yielding the exact global bound sum objective <= 45*alpha. The Sigma_v F_v = 0 identity row is consumed exactly like the Sigma g_v = 0 and Sigma h_v = 0 rows: the identity enters the search ONLY through the per-state interval endpoints and these two aggregate sign conditions; nothing else is trusted.

The cone basis is the full frozen basis plus the mixed row: every state carries all four interval pairs' endpoints, every certificate format carries all eight coefficients, and the dual discovery optimizes over all of them (gamma/delta/epsilon/zeta/eta/theta equal to 0 is an output, never a constraint — the frozen m2/m3/m4 rows are present and reused; the m4 result that the optimal certificates carried epsilon = zeta = 0 is DATA, not a hard-coded rule).

- [ ] **Step 1: Write failing synthetic-certificate tests**

Cover (eleven-field State construction; Synthetic states mirror the m4 module-level doctest shape):

    from mixed_deficiency_cone import State
    from search_mixed_cuts import (
        Certificate, PrimalWitness, exact_affine_bound,
        verify_certificate, verify_primal_witness,
    )

    # d=20, deficiency/excess_balance zero: the witness tests need a ZERO
    # balance column value on the single support state (sum b*x = 45*0 = 0),
    # mirroring the m4 synthetic-test convention.
    S = lambda g_lo, g_hi, h_lo, h_hi, f_lo, f_hi: State(
        20, 0, 0, 0, 0, g_lo, g_hi, h_lo, h_hi, f_lo, f_hi)

1. `test_rationalized_certificate_checks_every_state`: two states `S(-4, 2, -8, 1, -3, 5)`, `S(-2, 4, -6, 3, -7, 1)`; `exact_affine_bound(states, "total_deficiency")` verifies.
2. `test_tampered_alpha_is_rejected`: alpha minus one -> ValueError.
3. `test_tampered_alpha_plus_is_rejected`: alpha plus one -> ValueError (least-alpha exactness pin).
4. `test_negative_eta_is_rejected`: eta "-1" -> ValueError.
5. `test_negative_theta_is_rejected`: theta "-1" -> ValueError.
6. `test_negative_gamma_delta_epsilon_zeta_still_rejected`: one parametrized subTest per frozen sign gate.
7. `test_wrong_certificate_arity_is_rejected`: a dict or dataclass with fewer/more coefficient fields -> ValueError (schema-shape pin).
8. `test_exact_primal_rejection_witness_with_mixed_rows`: single state `S(-1, 1, -2, 2, -1, 1)`, witness {0: "45"} on "degree20_count" -> Fraction(45).
9. `test_primal_witness_violating_f_row_is_rejected`: `S(-1, 1, 1, 2, -1, 1)[f_lo -> 1]` i.e. state with `f_lo > 0`, witness {0: "45"} -> ValueError (analog of the m4 h-row pin).
10. `test_primal_witness_violating_f_hi_row_is_rejected`: state with `f_hi < 0` -> ValueError.
11. `test_wrong_objective_id_is_rejected`: "not_registered" -> ValueError.
12. `test_perturbed_f_window_moves_alpha`: the same cone with f_lo shifted by +1 changes the least alpha by the certificate's eta coefficient — pins that the F row is genuinely wired in (guards against an F-dead search implementation).

- [ ] **Step 2: Run tests and observe missing-interface failures**

```bash
cd /Users/jinleic/jinleic-workspace/math
./.venv/bin/python -m unittest r55/tests/test_search_mixed_cuts.py -v
```

Expected: import failure for `search_mixed_cuts`.

- [ ] **Step 3: Implement the exact search machinery**

Freeze before the first full run — routes and edges identical to m3/m4, frozen verbatim:

    total_deficiency: objective=state.deficiency,        accept upper_bound <= 315
    degree20_count:   objective=int(state.d == 20),      accept upper_bound < 1
    deficiency_ge8_count: objective=int(state.deficiency >= 8), accept upper_bound < 1
    required_local_family: UNAVAILABLE_NO_COVER_CERTIFICATE (no hash-pinned complete
                     gluing-cover manifest exists in this repository)

Dual discovery LP (SciPy/HiGHS, propose-only): variables (alpha, beta, gamma, delta, epsilon, zeta, eta, theta), bounds [(None,None), (None,None), (0,None)*6], objective [1, 0, ... 0], and one row per state

    [-1, -balance, -g_lo, g_hi, -h_lo, h_hi, -f_lo, f_hi] <= -objective(state)

with the row layout exactly in this coefficient order. Convert the returned (beta, gamma, delta, epsilon, zeta, eta, theta) to `Fraction` via the frozen denominator ladder (100, 1000, 10000, 100000, 1000000); derive alpha exactly as the max over states of `objective - beta*balance - gamma*g_lo + delta*g_hi - epsilon*h_lo + zeta*h_hi - eta*f_lo + theta*f_hi`. `verify_certificate` resolves the objective through the frozen registry, rejects unknown IDs, enforces the six sign gates, checks every state inequality in exact arithmetic, and returns 45*alpha.

Primal discovery for every route not accepted by an exact upper certificate — variables x_s >= 0, with the LP discovery layout of 2 equalities + 6 inequalities:

    maximize sum objective(s) * x_s
    subject to sum x_s = 45                                  (equality 1)
             sum excess_balance(s) * x_s = 0                 (equality 2)
             sum g_lo(s) * x_s <= 0                          (inequality 1)
             sum g_hi(s) * x_s >= 0                          (inequality 2)
             sum h_lo(s) * x_s <= 0                          (inequality 3)
             sum h_hi(s) * x_s >= 0                          (inequality 4)
             sum f_lo(s) * x_s <= 0                          (inequality 5)
             sum f_hi(s) * x_s >= 0                          (inequality 6)

Use HiGHS dual simplex only to identify a sparse basic support (weights above 1e-9). If exact reconstruction fails, form one deterministic candidate set from those indices and the 24 states with smallest (absolute_dual_slack, state_index) — SUPPORT_TOLERANCE = 1e-9, FALLBACK_SUPPORT = 24, MAX_SUPPORT = 8 (2 equalities + 6 inequalities span a basis) — and enumerate subsets of at most eight states with every active/inactive choice for the six inequality rows (64 choices), solving each basis by `Fraction` Gaussian elimination. Do not expand this set after seeing results: failure becomes MIXED_CERTIFICATION_UNRESOLVED (a solver status or rounded coefficient never proves a negative result).

`verify_primal_witness` independently requires: nonnegative canonical weights, total weight 45, zero excess balance, sum g_lo <= 0 <= sum g_hi, sum h_lo <= 0 <= sum h_hi, sum f_lo <= 0 <= sum f_hi, and recomputes the exact objective. A route is REJECTED_BY_EXACT_WITNESS only at value > 315 (total_deficiency) or value >= 1 (either count objective); otherwise CERTIFICATION_UNRESOLVED. No caller-provided lambda or stored boolean controls acceptance. Do not tune objective definitions, thresholds, the denominator ladder, or the support rules after seeing results.

Frozen constants: N_TARGET = 45, EXPECTED_STATES = 3215, ANALYSIS_FILE = "engstrom_identity.json".

- [ ] **Step 4: Run the cut search and persist exact decisions**

```bash
cd /Users/jinleic/jinleic-workspace/math
./.venv/bin/python r55/src/search_mixed_cuts.py \
  --analysis r55/data/engstrom_identity.json
```

Expected terminal contract is exactly one of:

    MIXED CUT SEARCH COMPLETE status=MIXED_ACCEPTED_CUT accepted=N (N >= 1)

    MIXED CUT SEARCH COMPLETE status=MIXED_NO_CUT_IN_FROZEN_BASIS accepted=0

    MIXED CUT SEARCH COMPLETE status=MIXED_CERTIFICATION_UNRESOLVED accepted=0

MIXED_NO_CUT_IN_FROZEN_BASIS requires an exact rejection witness for every instantiable route. MIXED_CERTIFICATION_UNRESOLVED exits nonzero and cannot trigger a scientific disposition. Every persisted search record contains its objective ID, exact upper certificate and bound when available, exact lower witness and value when available, derived route status, and no unchecked acceptance field. Schema 3's `searches` receives exactly four records; `disposition` advances exactly as the m4 artifact's did.

---

### Task 4: Disjoint checker and evidence-grade disposition

**Files:**
- Create: `math/r55/src/check_mixed_certificate.py`
- Create: `math/r55/tests/test_check_mixed_certificate.py`
- Modify: `math/r55/bench/spec.json`
- Modify: `math/PROGRESS.md`
- Modify: `math/r55/README.md`
- Modify: `math/r55/notes/structural_constraints_2026-08-15.md`

This task mirrors the landed `check_m4_certificate.py` discipline and must not begin before Task 3's artifact exists.

**Interfaces:**
- Consumes: frozen validation files, graph catalogs, frozen m3/m4 artifacts, and `data/engstrom_identity.json` (schema 3).
- Produces an independent pass/fail verdict. It must not import `subgraph_identities`, `m3_deficiency_cone`, `m4_subgraph_identities`, `m4_deficiency_cone`, `mixed_subgraph_identities`, `mixed_deficiency_cone`, or `search_mixed_cuts`. It may import only Python standard-library modules and `check_ramsey.parse_graph6_line`/`popcount`.

- [ ] **Step 1: Write checker corruption tests**

Create a temporary copy of the analysis artifact and independently test rejection of each of these twelve corruptions:

1. one changed source SHA-256;
2. one changed state endpoint (any of g/h/f lo or hi);
3. one reduced alpha coefficient;
4. one changed primal weight;
5. one false exact bound or route status;
6. one changed 656-graph replay count or n=49 replay field (`row_values`, `zero_corners`, `combos`, `extremal_graphs`);
7. one changed motif histogram count (any of the fifteen aggregates);
8. one changed outer-window endpoint (any aggregate);
9. one changed trust-root statement (any of the eight);
10. one deleted or added state record;
11. one swapped f-window pair (f_lo/f_hi exchanged on a state whose interval is not a point);
12. one changed schema_version (3 -> 2) or campaign id.

The untouched artifact must pass. Separately exhaust all labeled graphs through six vertices against the checker's own residual implementation, and show that the checker's kernel ground truths include the six corrected-math mutation pins — each wrong form fails on its named witness graph: (i) claw short form on the 4-vertex diamond (sum C(d,3) = 2, claw = 0); (ii) claw's missing -2*diamond silently broadens the diamond=792 X1/X2 strata; (iii) K1,3(bar) = paw on Z = K3+K1 (gives 0, truth 1: complement(K3+K1) is the claw, k3k1(K3+K1) itself is 1); (iv) C4(bar) = C4 on Z = 2K2 (gives 0, truth 1: two_k2(2K2)=1); (v) pc3 dropped-c_uv form on Z = K3 (gives -3, truth 0); (vi) K2+2K1 without the -2*(2K2) subtraction on Z = 2K2 (gives 2, truth 1) and paw without -12*K4 on Z = K4 (gives 12, truth 0).

- [ ] **Step 2: Implement the independent checker**

The checker deliberately reimplements this trust kernel from scratch:

- graph complement, induced subgraph, edge, wedge, triangle, induced-P3, independent-triple counts (its own naive + fast forms, cross-validated n <= 6);
- all eleven raw 4-vertex induced counters and the three pc3 closed forms, with the corrected universal identities (claw with -2*diamond -4*K4; K2+2K1 with the -2*(2K2) subtraction);
- the complement-pair map (K13(bar) = K3+K1; C4(bar) = 2K2; pc3 line; diamond(bar) = K2+2K1; paw(bar) = P3+K1; K4(bar) = i4; P4 self-complementary);
- `mixed_vertex_row` and `mixed_residual` recomputed as zero on all 328 published graphs and their complements, and the n=49 type-combo replay numbers (four combos of the two 132-edge census graphs; sorted values [0, 144, 288, 432]; exactly one zero corner);
- degree-histogram enumeration and the m3 q-window arithmetic for g endpoints; the m4 row endpoints h recomputed from the relation system and motif windows (as check_m4_certificate does);
- per-state f_lo/f_hi recomputed from the artifact's own motif windows through its OWN affine-interval code and the frozen Task 2 expansion (coefficients 1890/2109/135/-2/4124/-12/-528/12/4/-86/4 exact part; 516/72/24/24/24/(564-24d)/32 X-side; (74+4d)/-12/-8/-8/-24 Z-side) — the checker re-derives the F endpoints; it does not trust the producer's state records;
- the R1–R5 relation system with the m4-corrected constants and the outer-window derivation for the degree-histogram aggregates of q, t, p3, pc3;
- exact per-class motif histograms re-streamed from all applicable published-complete files (all fifteen aggregates);
- the cone-stability projection (rebuilt states match the frozen m4 artifact's 3215 records on (d, a, b, g_lo, g_hi, h_lo, h_hi));
- the frozen objective registry, eight-coefficient exact dual certificates, exact primal witnesses against the REBUILT states, bounds, thresholds, and derived route/disposition statuses;
- source hashes (frozen m3/m4 artifacts included), catalog counts, frozen replay counts, and the v3 key sets.

Envelope-LP scoping (verbatim from the m4 checker's document contract, kept deliberately): for the missing classes the checker trusts the PRODUCER's envelope-LP windows under caps + schema validation only; exact re-derivation there is scoped to its docstring as m4's was — a future MIXED_ACCEPTED-tier disposition under any follow-on campaign requires the plan-letter R1–R5 LP recompute before the evidence label holds.

Stage ordering (fail-fast, first mismatch prints the first exact counterexample and exits 1): schema/key-sets -> trust roots -> hashes and sizes -> input-set reconciliation -> catalog resweep (line counts, k4 == 0 per graph, fifteen aggregate histograms) -> outer-window audits -> state rebuild (g, h, then f) -> replay (656/656) -> n49 type-combo replay -> search records -> disposition.

Any producer count, bound, boolean, or disposition is recomputed, never trusted.

- [ ] **Step 3: Run all focused tests and the independent checker**

```bash
cd /Users/jinleic/jinleic-workspace/math
./.venv/bin/python -m unittest discover -s r55/tests -p 'test_*mixed*.py' -v
./.venv/bin/python r55/src/check_mixed_certificate.py \
  r55/data/engstrom_identity.json
```

Expected: focused tests `OK`; checker prints `MIXED EVIDENCE VERIFIED: <disposition>`.

- [ ] **Step 4: Apply the frozen scientific gate**

If at least one exact upper certificate meets a quantified success route:

- set `next_research_campaign.status` to `MIXED_ACCEPTED_PENDING_NOVELTY_AUDIT`;
- record the theorem statement with exact coefficients and hypotheses in the existing structural note;
- conduct the primary-source novelty audit (Engstrom arXiv:1002.4304 incl. Theorem 2.18/Corollary 2.19 and Conjectures 1.1/1.2; MR97 Conjecture 1, Lemma 2.3, Theorem 2.2/3.2 and section 5; Angeltveit-McKay arXiv:2409.15709) before using "new" or "publishable";
- do not claim R(5,5) <= 45.

If every instantiable route is rejected by an exact primal witness:

- set status to `MIXED_NO_CUT_IN_FROZEN_BASIS`;
- record each exact witness and state only that the frozen continuous basis {balance, g, h, F} cannot prove an accepted route;
- state explicitly that the campaign's first_exact_basis (MR97 m=2,3,4 plus the Engstrom mixed identity) is now exhausted, fulfilling stop condition 2 of `next_research_campaign`;
- preserve all artifacts and hashes;
- record the terminal fallback order for residual tiers — VeriPB-certified search, then srg(45,22,10,11) complementary SAT — as the decision point for the next separately reviewed campaign.

If any route is CERTIFICATION_UNRESOLVED:

- set status to `MIXED_CERTIFICATION_UNRESOLVED`;
- preserve the verified upper/lower bounds without a positive or negative claim;
- stop and record the exact unresolved gap and basis-reconstruction failure.

- [ ] **Step 5: Final integrity verification**

Parse the v3 artifact, verify every recorded SHA-256, rerun the checker, confirm the frozen m3/m4 artifacts, benchmark hashes, and candidate artifacts are unchanged, and rerun the m3, m4, and mixed focused suites. Report the mathematical result first, then files and verification evidence.

---

## Plan Self-Review

- **Spec coverage:** Implements the Engström mixed-identity tier of `next_research_campaign` under its `first_exact_basis` (the fourth listed identity) and its architecture slices: identity derivation and replay (Task 1, universal residual + 656-graph + n=49 type-combo ground truths), catalog motif vectors (Task 2, fifteen hash-pinned exact aggregates with a fail-closed streaming budget), the n=45 exact feasibility cone (Task 2, eleven-field State with the F interval layered on the frozen m2/m3/m4 rows), exact certificate search over the enlarged basis (Task 3, eight rational coefficients, 2-equality + 6-inequality discovery layout), and theorem extraction with a disjoint checker (Task 4). Stop-condition bookkeeping is explicit: a negative disposition here discharges stop condition 2 of the slot and records that fact; no pivot to search-engine optimization, SDP, or raw catalog scaling is introduced.
- **Placeholder scan:** No unresolved implementation placeholders, deferred code bodies, or production stubs remain. Implementation-time constants are confined to (a) the naive reference counters local to test files, (b) the mutation-selector coefficient tables in test files, and (c) internal factorizations of the closed forms — each pinned by an observable test contract, none a stub.
- **Type consistency (v3-typed):** `State` (eleven int fields; frozen dataclass, ordered), `Certificate` (objective_id + eight canonical reduced-fraction strings), `PrimalWitness` (objective_id + sorted index/weight pairs), the fifteen histogram keys, artifact path `data/engstrom_identity.json` at `schema_version: 3`, campaign id `higher_order_identity_positive_deficiency_mixed`, `EXPECTED_STATES = 3215`, `MAX_SUPPORT = 8`, `FALLBACK_SUPPORT = 24`, `SUPPORT_TOLERANCE = 1e-9`, and the CLI terminal contracts are consistent across all four tasks. The frozen m3/m4 artifacts are consumed read-only and never mutated; the m4 State's nine-field shape is never widened in place.
- **Scope:** No m=5 (void), no K>=5 Engström work (open conjectures), no SDP, no gluing engine, no VeriPB or SAT implementation, no benchmark or raw-generation work enters this plan; the fallback order is recorded, not executed.
- **Trust discipline:** The identity is re-derived in the module docstring ledger (five anchor lemmas) AND mechanically verified exhaustively n <= 6; the catalogs keep the four data trust roots of m3; the identity and the window tables are roots 7–8; AI/SciPy stay propose-only outside the trust base; every accepted number is rationalized and re-verified exactly, twice (producer + disjoint checker) with the frozen m4-artifact projection as a third anchor.

**Novelty table** (from the verified research inputs):

| Candidate deliverable | Closest precedent | Why it would be new |
|---|---|---|
| Exact mixed-identity certificate cone (T-family + Engström K=4 row) over complete published catalogs with exact rational certificates and a disjoint checker | MR97 §5, LP(s,t,n) — hand-derived, uncertified | first machine-certified exact basis containing a NON-separable subgraph-counting row; not in print |
| Engström K=4 identity inside an exact LP certificate at n=45 | Engström 2010 (identity only, no computational certificates); MR97 Conjecture 1 (unproven at the time) | the identity has never before been used in any computational certificate of any Ramsey claim |
| An accepted n=45 degree/deficiency restriction | none structural beyond the degree window {20..24} | first published n=45 structural restriction; directly shrinks the pointed-gluing lists behind R(5,5) <= 46 |
| Decisive basis-level negative with exact primal witnesses | none for this basis | publishable methodology evidence closing the ENTIRE recorded exact basis (m=2,3,4 + mixed) at n=45; the rank-deficit measurement (which rows the cone can and cannot see) is itself new evidence |

**If a cut accepts, the theorem is** (route-dependent, exact quantifiers): *For every graph G with |V(G)| = 45, omega(G) <= 4, and alpha(G) <= 4, it holds that* — `total_deficiency`: sum_{v} s(v) <= 45*alpha <= 315; `degree20_count`: min_v d_v >= 21; `deficiency_ge8_count`: s(v) <= 7 for every v — *conditional on the named catalog-completeness trust roots and the Engström K=4 identity, with the recorded rational certificate (alpha, beta, gamma, delta, epsilon, zeta, eta, theta) verified by the disjoint checker.* None of these claims R(5,5) <= 45.

**Research-input corrections applied during plan drafting** (each machine-verified in pre-plan session by Main or derivation-checked here; the plan's tests enforce them):

1. Complement map, two flipped lines vs the first brief (Main machine-verified 2026-08-18, 500 random labeled graphs per line, extraction note §5): K1,3(bar Z) = s(K3+K1, Z) — NOT paw (claw complement is triangle + isolate); s(C4, bar Z) = s(2K2, Z) — C4 is NOT self-complementary (its complement on four vertices is 2K2); only P4 is self-complementary. The two wrong candidate forms were themselves rejected by the same battery before landing.
2. pc3 (s(P3, bar Z)) requires the +c_uv term: pc3(Z) = Sum_{uv in E(Z)} (nu - d_u - d_v + c_uv); the dropped-c_uv form evaluates to -3 on K3 alone (witness Z = K3: truth 0). Two further provably equivalent closed forms are pinned as cross-checks: form B = e(nu-2) - 2*p3 - 3*t; form C = Sum_v C(nu-1-d_v,2) - 3*i3.
3. The claw closed form needs BOTH correction terms (this drafter's derivation, confirmed by Main on 600 random graphs n in 4..8 plus the pinned 4-vertex diamond witness): claw = Sum_v C(d_v,3) - paw - 2*diamond - 4*K4. The extraction note's original section-3 line carried only -4*K4; the relayed short form carried neither; both are captured as detected mutations in Task 1 with the diamond witness (degrees (3,3,2,2): Sum = 2, claw = 0).
4. K2+2K1 needs the -2*(2K2) subtraction: the bare closed-neighborhood binomial sum overcounts each induced 2K2 set twice (witness Z = 2K2: bare sum = 2, true induced K2+2K1 = 0); pinned in Task 1.
5. paw needs the universal -12*K4 term (witness Z = K4: triangle-degree sum = 12, paw = 0); on the swept Ramsey catalogs the k4 == 0 assertion fires first and then only does the term drop.
6. The n=49 Engström type-combo values {0, 144, 288, 432} (unique zero corner) are a kernel fixture only: they do NOT reproduce the MR97 Thm 3.2 contradiction (all values non-negative); the mixed row is a candidate CUT row at n=45, not a label/replay row.

## OPEN RISK: projection collapse (exit criteria)

The librarian's caveat is restated as an explicit open risk and bound to a concrete exit criterion. At the level of FORMAL identities the mixed row is disjoint-support and non-separable (MR97 Lemma 2.3 exhausts the separable identities of degree <= 6), but after projection onto the frozen (5,5,45) domain — d in {20..24}, K4-free N side, the R1–R5 relation system on the 4-motifs, and the frozen continuous relaxation — coordinate collapse is NOT ruled out: whether the mixed row increases the certificate cone rank is a COMPUTATION, not a formal corollary. The in-repo precedent is concrete and unfavorable: the m=4 tier added a genuinely new row under the identical discipline and the exact search still returned certificates with epsilon = zeta = 0 (zero strength gain); envelope sharpening, exact g-signs, census tightening, and parity refinements moved only route 2 and only modestly (notes/cone_lp_robustness_2026-08-17.md). It is entirely possible the F interval collapses into the shadow of the frozen rows and every optimal certificate returns eta = theta = 0. Exit criteria, both evidence-grade and both recorded: (a) NEGATIVE — the tier closes, the campaign's recorded first_exact_basis is exhausted (stop condition fulfilled), and the residual recorded fallback order (VeriPB-certified search, then srg(45,22,10,11) complementary SAT) becomes the decision point for the next separately approved design; (b) POSITIVE — MIXED_ACCEPTED_PENDING_NOVELTY_AUDIT gates on the primary-source novelty audit listed in Task 4, Step 4. Neither outcome may be relabeled: an eta = theta = 0 outcome on all accepted certificates is the specific, honest-negative shape the terminal record must name.

**Open questions for the main agent:**

- (a) State-domain cross-check strength: the plan RE-ENUMERATES the (d,a,b) domain and hard-asserts equality of the rebuilt 3,215-record projection (g/h endpoints and interval sources) with the frozen m4 artifact, rather than importing the m4 state records and only adding f endpoints. Recommendation: keep the assert (this plan); import-only would silently inherit any future m4-cone changes.
- (b) Envelope-audit scope for the eleven raw coordinates: compute LP outer windows for ALL classes (missing authoritative, present audit-only — this plan, mirroring the m4 choice) versus missing classes only. Recommendation: all classes; the audit-only recomputations double as relation-system validation on real data and cost the same basis enumerations.
- (c) Sweep width: the plan sweeps all eleven raw 4-motifs per graph (eleven histograms) rather than only the seven the F row's X side names, because the four complement-mapped motifs (two_k2, k2_2k1, p3_k1, plus pc3) and the per-graph k4 == 0 assertion are needed for Y-side windows anyway, and the P4-by-complement closed form then requires the remaining raw coordinates. Net effect: none of the fifteen aggregates could actually be skipped; the sweep is already minimal at this granularity. Flagged for confirmation, not re-opened.
- (d) The n=49 combo fixture ordering [0, 144, 288, 432] is pinned as the SORTED multiset of the four (X_i, Z_j) rows; the per-assignment mapping (which of the four type combos yields the unique zero corner) is intentionally NOT load-bearing for the cone — record it in the v3 artifact as supplementary evidence only if cheap.

**Fallback.** Recorded, not executed, under separate approval, in the terminal-record order for residual tiers after a negative disposition: (1) VeriPB certificates for the existing gluing runs; (2) complementary SAT on srg(45,22,10,11)-type configurations. This plan's tier occupies the middle slot of the m3/m4-era fallback list and is now the executed one; on MIXED_NO_CUT_IN_FROZEN_BASIS the terminal record states that the Engström tier is discharged and the basis-level stop condition of `next_research_campaign` is fulfilled.
