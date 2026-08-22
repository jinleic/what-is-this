# The Q=3, no-concurrency branch of K(14) is analytically CLOSED: no 54-family

**Verdict: OBSTRUCTION (proved).** *No arrangement of 14 distinct lines with
exactly three parallel pairs and no three concurrent lines carries a family of
54 pairwise interior-disjoint triangles.* Combined with the exact-verified
snapped Bader certificate (53 empty disjoint triangles, exactly Q=3), the value
of this branch of the K(14) case tree is **exactly 53**.

Everything below is machine-checked with exact rational arithmetic;
reproduction commands and artifact paths are in §9.

---

## 1. Setup and the equality regime

**Capacity theorem** (report.md §7): for lines with $Q$ parallel pairs and
points of multiplicity $k_p \ge 3$, any family $S$ of triangles with pairwise
disjoint open interiors satisfies

$$3|S| + C \;\le\; n(n-2) - 2Q - \sum_p k_p (k_p-4), \qquad C = \sum_T c_T,$$

with $c_T$ = number of non-supporting lines meeting the open interior of $T$,
and the per-line refinement $\sigma_\ell + \gamma_\ell \le m_\ell + 2a_\ell$.

At $n=14$, $Q=3$, **no multiplicities** ($k_p \ge 3$ absent), any 54-family
forces
$$3\cdot 54 + C = 162 + C \;\le\; 14\cdot 12 - 6 = 162 \quad\Longrightarrow\quad C = 0,$$
and **equality in every per-line inequality**: $\gamma_\ell = 0$ and
$\sigma_\ell = m_\ell$ for every line $\ell$, where $m_\ell$ = number of
*bounded elementary gaps* on $\ell$ ($11$ for the six paired lines, $12$ for
the eight free lines; total $6\cdot 11 + 8\cdot 12 = 162$).

So the **equality hypothesis** — the running assumption of everything below —
is:

> **(SAT)** every bounded elementary gap of $\mathcal A$ is the image of
> exactly one selected side.

**Lemma E0 (cells).** Under (SAT): every selected side IS one elementary gap,
and every selected triangle is a cell (empty zigzag triangle).

*Proof.* If side $[P_iP_j]$ on $\ell$ contained an interior vertex $P_k$,
the line through $P_k$ is transversal to $\ell$ and hence meets the open
interior of the triangle, contradicting $C=0$. A triangle whose sides are
gap-segments and whose interior meets no line is exactly a face. $\square$

**Lemma E1 (used gap = unique triangular side).** Under (SAT): every gap $e$
has **exactly one** adjacent triangular cell (its user). The other adjacent
cell is never a triangle: with $\ell$ placed horizontally and $c(P_i)$ = cot
of the $i$-th crossing line relative to $\ell$, the cell above
$e_t = [P_tP_{t+1}]$ is a triangle iff $c_t > c_{t+1}$ **and** the apex
triangle $(\ell, f_t, f_{t+1})$ is uncrossed; below iff $c_t < c_{t+1} \ \&$
uncrossed. The two conditions are mutually exclusive ($c_t \neq c_{t+1}$ in
the simple regime; equality happens only for a *flat gap* — see L5).

*Verified exactly*: on every gap of every analyzed arrangement (the cot/side
convention `assert (kind=='D') == (apex above)` never fired).
∎ (gap/apex convention certified on 162+162 Bader gaps and all small cases)

**Divisibility.** Under (SAT), used gaps $= 3 \cdot \#\text{cells}$
(each cell uses 3 gaps, each gap has exactly one user) so
**escapes $\equiv 0 \pmod 3$ at $(n,Q)=(14,3)$**: the smallest nonzero number
of escapes is 3.

## 2. The vertex analysis (L2, L3): alternation and extremity

**Lemma L2 (zigzag).** Under (SAT), at every *interior* vertex of every line
the cot sequence alternates strictly: $c_{t-1} < c_t > c_{t+1}$ or
$c_{t-1} > c_t < c_{t+1}$.

*Proof.* At $P = \ell \times f$ the four gaps at $P$ ($e^-, e^+$ on $\ell$;
$g^-, g^+$ on $f$, cyclically) bound four quadrant cells $x_1..x_4$. Each
bounded gap has exactly one triangular adjacent cell (L1); a triangular
quadrant cell covers exactly the 2 gaps of its quadrant. The cyclic XOR
$t_1+t_2 = t_2+t_3 = t_3+t_4 = t_4+t_1 = 1$ has only the solutions
$(1,0,1,0)$, $(0,1,0,1)$: *two opposite cells are triangles*. Quadrant "NW"
is triangular iff $c_{t-1} > c_t$, etc.; the two solutions read out as
"$c$: ascent-then-descent" or "descent-then-ascent". $\square$

**Lemma L3 (no mixed extreme; exactly n double-extreme vertices).**
Call a vertex $P = \ell\times f$ **extreme** on $\ell$ if it is the first or
last crossing on $\ell$, and **double-extreme** if extreme on both. Under
(SAT):

* **(L3a)** no vertex is extreme on exactly one of its lines;
* **(L3b)** the number $c$ of double-extreme vertices equals $n = 14$.

*Proof.* (a) Suppose $P$ is rightmost on $f$, interior on $\ell$. The two
cells adjacent to the empty $f$-ray at $P$ are unbounded, hence not
triangles. Both $\ell$-gaps at $P$ are bounded, so by L1 each has a triangular
cell, forced to be the *other* cells $X_0, X_1$ adjacent to them — and both of
*those* are adjacent to the bounded $f$-gap: both triangle-shaped, contradicting
L1's one-sidedness. (b) Every line has exactly two extreme vertices; extreme
incidences count $2n = 2c + b$ with $b$ = number of mixed vertices; (a) gives
$b = 0$, so $c = n$. $\square$

**Lemma L4 (direction adjacency).** Under (SAT), at every double-extreme vertex
$P = \ell \times f$ the directions of $\ell$ and $f$ are **adjacent in the
circular order of slope classes**: no line's direction lies strictly inside
the cone of the two empty rays.

*Proof.* The wedge $W$ between the two empty rays is a convex open cone
$P + \mathrm{cone}(u,v)$, angle $< \pi$. A line $g$ with direction strictly
inside $\mathrm{cone}(u,v)$ crosses both support lines $\ell, f$; its
crossings must lie on the *non-empty* rays (else $P$ isn't extreme), but a
line with direction inside the cone through two points on the opposite rays
never enters the cone's open side: an affine two-line computation
($g$ in cone coordinates passes through $(-\alpha,0),(0,-\beta)$, hence misses
the open quadrant) — contradiction. $\square$

**Lemma L5 (flat gaps).** A *flat gap* — a gap whose two endpoints are the
crossings of $\ell$ with the two members $a \parallel b$ of one parallel pair
and with NO other crossing in between — has both adjacent cells unbounded
(half-strips), so it can never be used. Under (SAT) flat gaps do not occur: on
every $\ell$, some third line crosses strictly between $\ell \cap a$ and
$\ell \cap b$, for every parallel pair $(a,b)$. (This is the correct
"the injection's tokens are unused" refinement of the task's candidate lemma:
in the simple regime D$\wedge$A is *automatically impossible* — a gap is a
strict descent XOR a strict ascent XOR flat; nothing is "both".)

**Lemma L6 (status of "no gap is both D and A").** With no concurrency, no
elementary gap is both D and A: at each endpoint of $e_t$ a *single* line
crosses $\ell$, so the cot comparison is total; a gap is a strict descent, a
strict ascent, or flat. The candidate consequence is TRUE but vacuous in this
branch (verified on Bader: count of both-D-and-A gaps $= 0$). Bader's
certificate is NOT in the equality regime (159 used, 3 escapes, $b=2$), so no
lemma with hypothesis (SAT) must — and does — cover it.

## 3. THE OBSTRUCTION: the extremity multigraph

**Theorem M.** No arrangement of 14 lines with exactly three parallel pairs,
no three concurrent, satisfies (SAT).

*Proof.* Assume (SAT). By L3 there are $14$ double-extreme vertices; by L4
each sits at the crossing of two **cyclically adjacent** slope classes. Let
the slope classes be $C_1,\dots,C_{11}$ in circular order ($11 = 14 - 3$) with
sizes $s_i \in \{1,2\}$, three of them $2$ — the alternative layout, one class
of size 3 ("triad"), is treated separately below. Let
$$m_i \;=\; \#\{\text{double-extreme vertices at } C_i \times C_{i+1}\},\qquad i \in \mathbb Z_{11}.$$

Every line of class $C_i$ contributes its two extreme vertices, which are
double-extreme (L3a) and adjacent-class (L4), so

$$\text{(M1)}\quad m_{i-1} + m_i = 2\, s_i, \qquad
\text{(M2)}\quad 0 \le m_i \le s_i\, s_{i+1}, \qquad m_i \in \mathbb Z .$$

(M2) is the crucial *multiplicity bound*: classes $C_i, C_{i+1}$ have
$s_i s_{i+1}$ mutual crossings in total, and distinct double-extreme vertices
are distinct crossings (each one uses its own pair of lines).

Now propagate (M1) along a run of $r$ consecutive single classes between two
doubled classes. The singles' equations $m_{j-1} + m_j = 2$ force strict
alternation $m \in \{a, 2-a\}$. Since single-single edges satisfy
$m \le s_i s_{i+1} = 1$ (two single lines cross *once* — they cannot
contribute two double-extreme vertices), every run of singles of length
$\ge 3$ has all multiplicities $= 1$; checking the finitely many compositions
$r_1 + r_2 + r_3 = 8$ with the endpoint equations $m = 4 - (\text{run-end})$
at each doubled class against (M2) fails for every layout; the unique exact
rational solution of the (full-rank, odd-cycle) linear system is **infeasible
for all $\binom{11}{3} = 165$ patterns** (66 have a negative entry, 99 force
$m_i \ge 2$ onto a single-single edge; zero are integral-feasible).

For the **triad** layout — one class of size 3 plus **eleven** singles
(14 lines, $R=12$, an EVEN cycle) — (M1) requires the even-cycle consistency
$\sum_i (-1)^i\, 2 s_i = 0$; with one 3-class the sum equals $\pm 4 \neq 0$ for
every one of the 12 placements: **no solution, not even rational**. (An
earlier draft of this certificate used an 11-class/13-line stand-in for the
triad by mistake; the 12-class/14-line leg above is the certified object —
thanks to reviewer Main for the catch.)

By (M1)+(M2) infeasibility in both layouts, no saturated arrangement exists.
$\square$

**Corollary.** Any 54-family at $n=14$ with $Q=3$ and no 3 concurrent lines is
impossible. With the divisibility note, the next best is $\ge 3$ escapes, i.e.
$|S| \le 53$; Bader's snapped certificate attains $53$: **branch value = 53.**

**Consistency controls (must stay feasible):** the SAME system is feasible for
every known equality arrangement: pentagon/triangle/Blanc-9/rot-sym-15 give
$m = (1,1,\dots)$; the affinely regular hexagon (6 lines, 3 parallel pairs,
exactly this layout with $R=3$) gives $m = (2,2,2)$ — equality WITH three
parallel pairs EXISTS at $n=6$, so the theorem is genuinely $n=14$-specific.
(M-checked and adjacency-checked exactly.)

## 4. Exact verification on Bader's certificate

Certificate: `scratch/kobon/discoveries/n14_bader53_verified.json`,
`savchuk-stack-straightened`, 53 triangles listed as line triples, exact
rational coefficients; lines $\ell_i: y = a_i x + b_i$.

* **Stored coefficients** carry **two** exact parallel pairs:
  $(\ell_0,\ell_1)$, $(\ell_2,\ell_3)$; lines 6, 7 are parallel to
  $9\cdot 10^{-15}$ relative accuracy but not exactly. Stored stats:
  $Q=2$, simple, 164 gaps, **159 used, 5 escapes**.
* **Snapped variant** ($\ell_7$ slope $:= \ell_6$ slope, intercept kept):
  $Q=3$, simple: re-verified the whole family exactly:
  all 53 triples non-parallel, all 53 triangles **empty** (no other line
  crosses any open interior, all $\binom{53}{2} = 1378$ pairs open-disjoint
  by exact Sutherland–Hodgman clipping), every side an elementary gap, the
  side-to-gap map injective: **159 gaps used, 162 total, escapes exactly 3.**
* Both-D-and-A gaps: **0** (L6 verified); flat gaps: **0** (L5 consistency).
* The three escapes are *nonempty zigzag triangles* (L1's shape exists, the
  apex-triangle is crossed):

  | line | slot | kind | gap endpoints (crossing classes) | cutting lines |
  |---|---|---|---|---|
  | 7 | 9 | D | crossings with $\ell_{11},\ell_{9}$ | 1, 3, 8 |
  | 10 | 11 | D | crossings with $\ell_{12},\ell_{11}$ | all 11 others |
  | 11 | 11 | A | crossings with $\ell_{7},\ell_{10}$ | 0,1,2,3,4,12,13 |

  They cluster at Bader's **two mixed vertices** (7,9) and (10,12):
  $b = 2$, $c = 13$, $2c + b = 28 = 2n$ — exactly the forbidden-extremity
  pattern of L3, failing the (SAT) hypothesis as it must.
* All 13 double-extremes of snapped-Bader satisfy direction-adjacency L4.

## 5. Validation ledger on known arrangements

| arrangement | n | Q | gaps | used | cells | esc | zz-viol | b | c | eq? |
|---|---|---|---|---|---|---|---|---|---|---|
| triangle | 3 | 0 | 3 | 3 | 1 | 0 | 0 | 0 | 3 | YES |
| pentagon side-lines | 5 | 0 | 15 | 15 | 5 | 0 | 0 | 0 | 5 | YES |
| affine hexagon side-lines | 6 | 3 | 18 | 18 | 6 | 0 | 0 | 0 | 6 | YES |
| heptagon side-lines | 7 | 0 | 35 | 21 | 7 | 14 | 14 | 0 | 7 | no |
| nonagon side-lines | 9 | 0 | 63 | 27 | 9 | 36 | 36 | 0 | 9 | no |
| **Blanc 21 (gallery, exact)** | 9 | 0 | 63 | 63 | 21 | 0 | 0 | 0 | 9 | **YES** |
| record 32 (gallery, exact) | 11 | 0 | 99 | 96 | 32 | 3 | 6 | 6 | 8 | no |
| **rot-sym 65 (gallery, exact)** | 15 | 0 | 195 | 195 | 65 | 0 | 0 | 0 | 15 | **YES** |
| **Bader snapped 53 (exact)** | 14 | 3 | 162 | 159 | 53 | 3 | 4 | 2 | 13 | no |
| circle-tangent n = 14, Q=3, ALL $\binom{11}{3}=165$ doubling patterns | 14 | 3 | 162 | 42 | 14 | 120 | 120 | 8 | 10 | no |
| parabola-tangent n = 4..17 (debug family) | — | 0 | n(n−2) | 3(n−2) | n−2 | ... | ... | ... | 3 | no |

Observations the table makes precisely checkable:
* The equality regime occurs iff L2/L3/L4 all hold (`zz=0, b=0, c=n`).
* Circle-tangent (polygon/Cramer) types with 3 antipodal pair classes give
  **exactly 14 cells = n** at n=14 Q=3 for **all** 165 doubling patterns –
  the polygon world is far from the capacity (needs 162 = 3·54).
* The records that DO hit equality (n = 5, 6, 9, 15) are non-generic /
  high-symmetry types; at n=14 no pattern within reach of polygon methods works.

## 6. Support experiments (NOT part of the proof)

* Random exact-rational Q=3 arrangements (1401 samples, slopes/intercepts
  small rationals, 3 exact pairs): **zero** theorem violations; escape
  histogram $\in [45, 111]$; random landscape is far from the capacity.
* Bader-neighborhood hill climb (1200 exact local moves: intercepts,
  pair-slopes, single slopes): **3 stays locally optimal**.

## 7. Precise hypothesis ledger (task requirement)

| lemma | statement | hypothesis | status |
|---|---|---|---|
| E0/E1 | sides = gaps, triangles = cells, unique user per gap | (SAT)=$C{=}0$ + tight | PROVED here |
| L2 | zigzag cot at every interior vertex | (SAT) | PROVED; exactly verified on all cases above |
| L3 | no mixed-extreme vertex; $c = n$ | (SAT) | PROVED; $\checkmark$ hexagon/n9/n15 |
| L4 | direction-adjacency at double-extremes | simple only (wedge) | PROVED; $\checkmark$ hexagon/n9/n15/Bader |
| L5 | flat gaps unusable | simple only | PROVED |
| L6 | no gap is both D and A | simple only | true but vacuous-here; $\checkmark$ Bader: 0 |
| M | (SAT) impossible at (14, Q=3, simple) | n=14, Q=3, simple | PROVED; exhaustive check |

Bader is safe under every lemma: it violates (SAT) (3 escapes, $b=2$).

## 8. Consequences for the K(14) case tree

* `q3paired` cube: **closed analytically** (paired layout of M).
* `q3triad` cube: **closed analytically** (triad layout of M).
* Branch value: **53 exactly** (≤ by Theorem M, ≥ by the snapped Bader
  certificate, no SAT cubes needed either direction in this branch).
* The same M-machinery generalizes: at $n=18$, $Q=3$ the direction-cycle has
  $r=16$ (even) and the necessary alternating-sum consistency
  $\sum_i (-1)^i s_i = 0$ cannot hold for any placement of the three pairs
  (parity: the sum is an odd multiple of 2) — the q3 branch of n=18/T=94
  dies the same way (one-line computation).
* Remaining legs: Q ≤ 2 and Q + multipoints (≠ this branch).

## 9. Artifacts and reproduction

Directory `scratch/kobon/n14/escape/`:

* `escape_lib.py` — exact rational arrangement engine (vertices, gaps,
  D/A classification, apex/side conventions, emptiness, vertex statistics,
  Sutherland–Hodgman disjointness). All assertions passed everywhere.
* `escape_bader_check.py` → `bader_stored_report.json` + `bader_snapped_report.json`:
  the entire §4 ledger (family re-verified, escapes classified).
* `escape_small_cases.py` → `small_cases2.json` (circle-tangent series and
  even-m antipodal families).
* `escape_n14_q3_search.py` → `n14_q3_pattern_search.json` (165/165 patterns,
  escapes = 120 always).
* `escape_multigraph.py` → `multigraph_feasibility.json`: the M-system with
  exact Gaussian elimination: 0/165 feasible; controls feasible.
* `escape_theorem_verify.py`: V1 direction-adjacency on all five real
  arrangements (0 violations), V2 infeasibility both layouts, V3 controls.
* `escape_random_search.py` → `random_stress.json` (1401 samples, 0 violations).
* Bader-neighborhood probe: `best neighborhood escapes: 3`.

Commands (venv python =
`/Users/jinleic/jinleic-workspace/scratch/kobon-audit/venv/bin/python`):

```
cd scratch/kobon/n14/escape
$PY escape_bader_check.py
$PY escape_small_cases.py
$PY escape_n14_q3_search.py
$PY escape_multigraph.py
$PY escape_theorem_verify.py
$PY escape_random_search.py
```

*Prepared by agent EscapeObstruction, 2026-08-18; the capacity theorem it uses
is the 2026-08-17 crossing-refined theorem, VALID + novel; the obstruction
proof here is new and all verification is exact.*
