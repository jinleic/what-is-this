# The Kobon Triangle Problem: A Certified Exact Decision at n = 10

Computational-proof report — opened 2026-08-13, certified complete
2026-08-15. Code: `engine.py`; instances: `kobon_n10_t26_proof_*.cnf`;
certification audit: `AUDIT.md`.

**Theorem.** $K_{\rm gen}(10) = 25$: among all arrangements of 10 distinct
lines in the real plane — parallels, multiple intersection points, and
crossed triangles allowed — the maximum number of pairwise interior-disjoint
triangles with all three sides on the lines is exactly 25.

## 1. Problem

The Kobon triangle problem (Fujimura 1978; on Wikipedia's list of unsolved
problems in mathematics) asks for the maximum number of non-overlapping
triangles whose sides lie on $n$ straight lines. This report writes
$K_{\rm gen}(n)$ for the broad convention it decides: selected triangles have
pairwise disjoint interiors, other lines may cross a selected triangle, and
parallel lines and multiple intersection points are allowed.

**Target.** Decide whether $K_{\rm gen}(10)=25$. The lower bound is an
explicit integer line arrangement. The upper-bound campaign uses an exhaustive
family of necessary-clause CNFs with independently checked DRAT proofs. The
all-$n$ Kobon problem remains outside this finite decision.

## 2. State of the art — and a literature gap

Two distinct quantities appear in the literature:

* $as_3(n)$ — max **triangular faces** in a **simple** arrangement (every pair
  crosses, no three concurrent): BBL (arXiv:0706.0723), Blanc (arXiv:0801.2845).
  Sharp even-$n$ results: $as_3(8)=14$, $as_3(10)=25$, $as_3(12)=37$, $as_3(14)=53$.
* $K_{\rm gen}(n)$ — the broad Kobon quantity fixed in §1. The public
  constructions $K(8)\ge15>as_3(8)$ and $K(12)\ge38>as_3(12)$ show that
  arbitrary Kobon triangles can beat the simple triangular-face count.

**Scope gap (source audit, 2026-08-13).** Public tables record 25 and 38 as
exact at $n=10,12$, but the proof chain located for those entries mixes
different conventions:

* BBL/Blanc prove statements about triangular **faces** of **simple**
  arrangements. BBL explicitly say their even-$n$ bound is not attained for
  12 and 14 pseudolines.
* The Clément–Bader draft develops its charging argument from
  pairwise-intersecting, nondegenerate/perfect configurations.
* Savchuk 2025 can represent parallels and multiple points in pseudoline
  tables, but gives no exact theorem for $n\in\{10,12,14\}$.
* Parpalak–Utkin 2026 concerns bounded triangular faces in general position.

Thus those sources do not by themselves prove the broad convention fixed in
§1. When the Clément–Bader/Tamura local charging steps are extended to that
convention, exact rational counterexamples appear (§6). This report does not
claim the public numerical entries are wrong; it supplies a self-contained
certificate that degeneracies and crossed triangles cannot improve the
$n=10$ value.

For orientation, the located construction/bound record was:

| $n$ | construction | broad-convention ceiling asserted in public tables/drafts |
|----:|:----:|:----:|
| 10 | 25 | 25 (public table), 26 (Clément–Bader draft formula) |
| 12 | 38 | 38 (public table), 39 (Clément–Bader draft formula) |
| 14 | 53 | 54 (public table), 55 (Clément–Bader draft formula) |

The exact question decided here is whether ten lines in the broad convention
can form 26 pairwise interior-disjoint triangles. UNSAT proves they cannot;
the explicit 25-triangle configuration then determines the maximum.

## 3. Method

SAT decision procedure over a **combinatorial relaxation** of line configurations.

### 3.1 Normalization (WLOG)

Rotate so no line is vertical and all crossing x-coordinates along each line are
distinct per point; index lines by slope (ties = parallel classes, contiguous by
construction). Any configuration attains this after rotation + relabeling.

### 3.2 Variables

$P(i,j)$ (lines cross), $C(i,j,k)$ (concurrent triple), $X(r;i,j)$ (crossing of $i$
strictly precedes crossing of $j$ along line $r$), $U(i,j)$ (of two parallel lines,
lower index lies above), $S(t)$ (triple $t$ selected as a Kobon triangle), plus
defined auxiliaries $SA/SB$ (vertex strictly above/below a line), $B$ (betweenness
along a line), $IN$ (vertex strictly inside a triangle).

### 3.3 Machine-verified lemmas (exact rational arithmetic)

* **L1 (middle crossing).** For slopes $m_a<m_b<m_c$ pairwise crossing:
  $x_{ac}$ lies strictly between $x_{ab}$ and $x_{bc}$.
  Proof: $f=\ell_b-\ell_a$ increasing (zero at $x_{ab}$), $g=\ell_b-\ell_c$
  decreasing (zero at $x_{bc}$), $h=f-g$ zero at $x_{ac}$; sign analysis.
  Verified on 11,810 triples embedded in degenerate configurations: 0 failures.
* **L2 ($\sigma$ three-way equality).** One bit per triple:
  $X(a;b,c)=X(b;a,c)=X(c;a,b)$.
* **L3 (sidedness table).** $\sigma_{abc}=+ \iff P_{bc}$ above $a$ $\iff P_{ab}$
  above $c$ $\iff P_{ac}$ **below** $b$. (Ties $\to$ both order literals false,
  encoding "on the line" for free.)
* **L4 (parallel linkage).** For $u\parallel w$ and $v$ crossing both: $v$ meets
  the lower of $\{u,w\}$ first iff $m_v > m_u$.
* **L5 (disjointness rules R1–R4** $\equiv$ **exact interior-disjointness).**
  - R1: two selected triangles sharing side-line $\ell$ with interiors on the
    same side have side-intervals with disjoint interiors (opposite sides may
    share a segment — C–B Fig. 2).
  - R2: no vertex of one strictly inside the other.
  - R3: no transversal side×side crossing interior to both sides.
  - R4: a vertex lying (via concurrency) strictly inside another triangle's side
    forces the rest of that triangle weakly to the exterior side.
  Each rule has a local convex-geometry necessity proof. Equivalence with the
  exact separating-axis oracle verified on **167,557 random triangle pairs** in
  configurations with planted parallels/concurrencies: 0 necessity violations,
  0 sufficiency gaps.

### 3.4 Soundness architecture

Every clause is *necessary* (holds for every real configuration), so
**UNSAT at target $T$ proves $K(n)<T$** for arbitrary line configurations.
A SAT model is only a candidate: it is realized by fixed-slope LP straightening
(order constraints are linear in intercepts once slopes are fixed; concurrencies
are linear equalities solved exactly over $\mathbb{Q}$), then **verified with
exact rational arithmetic** — the resulting configuration is the lower-bound
certificate, independent of the model's correctness.

### 3.5 Proof accelerators (all necessary)

* **Exactly $T$ selections.** A family of more than $T$ disjoint triangles has
  a $T$-element subfamily, and deselection only disables clauses. The proof CNF
  therefore adds $\sum_t S(t)\le T$ to the base $\sum_t S(t)\ge T$ constraint.
* **Exact bounded-face penalty.** For an essential arrangement of distinct
  affine lines, the region count and $2n$ unbounded regions give
  \[
  F_b=1-n+\sum_p(k_p-1)
     =\binom{n-1}{2}-Q-\sum_p\binom{k_p-1}{2}.
  \]
  Every selected triangle contains a bounded arrangement face; disjoint
  interiors require distinct faces. Hence
  \[
  T\le\binom{n-1}{2}-Q-\sum_p\binom{k_p-1}{2}.
  \]
  The CNF represents each $k$-line point once by its three least-labelled
  lines and $u=k-3$ later incident lines, charging
  $1+2u+\binom u2=\binom{k-1}{2}$. An independent reviewer proved the formula
  including multiple vertices, crossed triangles, shared boundaries, and
  parallel classes, then audited every Tseitin clause and counter weight.
* **No-concurrency side capacity.** On line $\ell$, order its
  $m=n-1-q_\ell$ finite crossing points and write the crossing lines as
  $x=x_i+c_i y$. A triangle base $[i,j]$ is above $\ell$ exactly when
  $c_i>c_j$, and below exactly when $c_i<c_j$. Same-side bases have disjoint
  interiors; assigning above bases to descent gaps and below bases to ascent
  gaps injects all of them into the $m-1$ adjacent gaps. Thus
  $\sigma_\ell\le n-2-q_\ell$ and
  $3T\le n(n-2)-2Q$. At $(n,T)=(10,26)$ this forces $Q\le1$.
* **Exhaustive proof cubes.** The final family is: simple; one parallel pair
  with no concurrency; a normalized parallel pair with some concurrency; and,
  when every pair crosses, one concurrency triple in each of the eight
  $D_{10}$ orbits. Only $P/C$ units are fixed. No second orientation bit is
  fixed.

### 3.6 Validation ladder

The engine reproduces the public small-$n$ values, including cases where a
simple triangular-face count is too small:

| $n$ | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|---|
| $K(n)$ | 1 | 2 | 5 | 7 | 11 | 15 | 21 |
| engine | SAT@K, UNSAT@K+1 for all seven values | | | | | | |

End-to-end pipeline test at $n=8$: SAT model (2 concurrent triples) →
LP straightening → exact verification of 15 pairwise-disjoint triangles. ✔

## 4. Result

### 4.1 Exact lower-bound certificate

For line $i$, let $\ell_i:y=m_i x+b_i$ with

```text
(m_i,b_i) =
(-39,0), (-39,-9653), (-37,-7868), (-10,-589), (31,-2168),
(33,-2158), (51,-2878), (53,92), (58,0), (100,490).
```

Select these 25 triples of line indices:

```text
026 035 038 045 046 078 079 123 147 157 158 168 169
237 247 248 258 259 269 346 356 379 389 459 489
```

All ten affine lines are distinct; $\ell_0\parallel\ell_1$. Every listed
triple is nondegenerate. An exact `Fraction`/integer-homogeneous
separating-axis check verifies all $\binom{25}{2}=300$ triangle pairs have
disjoint interiors. `verify_n10_lower_bound()` returns `(True, "ok")`.
Therefore $K_{\rm gen}(10)\ge25$.

### 4.2 Upper-bound certificate — COMPLETE (2026-08-15)

The 2026-08-15 independent byte-level audit (`AUDIT.md`) found nine of the
eleven historical cube proofs structurally complete and two (`c0`,
`par_conc`) truncated when their original solver processes died. Resolution,
all independently machine-checked:

* **9 shipped proofs re-verified**: drat-trim backward checking returned
  `s VERIFIED` for `simple`, `par_noc`, `c1`–`c7`
  (≈23.4 CPU-h; `scratch/kobon-audit/results_verified.tsv`).
* **2 truncated cubes re-proven by verified case split**: each cube CNF was
  split on 3 crossing-order variables into 8 subcubes — each subcube file
  byte-verified as the parent CNF plus exactly 3 unit clauses, so the 8
  polarities form a tautological case split. kissat proved all 16 subcubes
  UNSAT (exit 20) and drat-trim verified all 16 fresh DRATs
  (`s VERIFIED`; ≈16.4 CPU-h; 44 GB of proofs;
  `scratch/kobon-audit/splits.tsv`).

Every cube of the exhaustive 11-cube cover therefore carries an
independently verified UNSAT certificate at target 26. With §4.1,
$K_{\rm gen}(10)=25$. $\blacksquare$

## 5. Notes from the attack

* **Certification pipeline validated**: $K(9)\le21$ UNSAT rerun with DRAT logging
  (351 MB proof, 120M resolution steps) and independently verified by `drat-trim`:
  `s VERIFIED` (444 s).
* **Exact auxiliaries pay**: making SA/SB/B reified (two-way) instead of
  one-directional cut the $n{=}9$@22 UNSAT time from 845 s to 352 s (kissat,
  same hardware, `kobon_n9_t22_exact.cnf`).
* **Cube-and-conquer**: the $n{=}10$@26 question splits soundly into
  simple / parallel / concurrency cases. The final **simple cube is UNSAT**:
  no simple 10-line arrangement admits 26 pairwise interior-disjoint triangles
  even when other lines may cross them — stronger than the face-count statement
  $as_3(10)=25$.
* **A cautionary lemma**: a natural Clément–Bader-style budget bound
  ($2\cdot\#\text{par-pairs} + \#\text{triple-pts} + \text{excess} \le$ slack,
  via "shared segment ⇒ multipoint endpoint") is **false** under the Kobon
  convention: machine search found 256/500 random configurations with segments
  double-covered from opposite sides and no multipoint — two triangles' side
  intervals may partially overlap, their vertices subdividing each other's sides
  (excess pays for sharing). Exact rational counterexamples later showed why
  those charging steps cannot simply be extended to the broad convention (§6).
  Budget-restricted cubes were **discarded as non-proof-relevant**; only
  assumption-free instances count.
* **Independent adversarial audit** (reviewer agent, full pass over A1–A7,
  R1–R4, auxiliaries): *no unsound clause found* — replicated necessity on
  1,296 exact-rational configurations (planted parallels/concurrencies, 0
  violations) and 203,673 triangle pairs against the separating-axis oracle
  (0 gaps either direction). Three real defects surfaced and fixed:
  1. *Symmetry-breaking justification*: the "x-flip" argument was false (the
     x-flip acts as identity on $\sigma(0,1,2)$ at $n{=}3$). The clause is
     sound via the $\pi$-rotation $(x,y)\mapsto(-x,-y)$ ($b\mapsto-b$), which
     preserves labels and flips **all** $\sigma$ bits — so **exactly one**
     $\sigma$ bit may be fixed WLOG. An explicit 4-line witness shows fixing a
     second bit is unsound. A cube family that added per-cube $X$-units on top
     of `add_symbreak` was consequently **killed before use**.
  2. *Reproducibility*: proof CNFs must be regenerable from source; `engine.py`
     now contains the full instance driver (`dump_instances`), and all 12
     base/cube CNFs regenerate **byte-identically** in a fresh environment.
  3. *Dihedral WLOG has a twist*: a plane rotation shifts slope-sorted labels
     cyclically ($i\mapsto i{+}k$) but XORs each line's order relation with
     $[d_r<0]$ and swaps SA/SB on flipped lines; naive lex-leader breaking over
     index permutations would be unsound. Concurrency variables $C$ permute
     *without* twist, so cube units on $C$-orbit representatives are safe.
* **Sound $D_{10}$ cover** for $n=10$@26 (each cube = base + necessary units;
  `add_symbreak` remains sound inside every cube since $b\mapsto-b$ is
  label-preserving):
  $\{$ *simple* (all $P$, no $C$); *par-noc* (exactly the normalized parallel
  pair $\ell_0\parallel\ell_1$, no $C$); *par-conc*
  ($\ell_0\parallel\ell_1$ and some $C$); *c₀..c₇* (all $P$ plus one of the
  eight dihedral-orbit representatives of concurrent triples) $\}$.

## 6. Why classical charging does not cover the broad convention

A reviewer attempted to extend the budget inequality
$3T \le n(n-2) - 2Q - \sum_p (k_p-2)^2$ to §1's convention. Exact rational
witnesses refute the required local steps when crossed triangles and
concurrencies are admitted. This does not refute a theorem stated only for
nondegenerate/simple configurations; it blocks using that argument as the
general upper-bound proof needed here:

1. **The extension of C–B Lemma 1 step (2) is false** ("a segment siding two triangles has a
   multiple point at an endpoint"): $n=4$, lines $y=0$, $y=x/5$, $y=x-1$,
   $y=2x-4$; triangles $\{0,1,2\}, \{0,2,3\}$ are interior-disjoint and share
   the whole elementary segment $(1,0)$–$(5/4,1/4)$ of $y=x-1$, both of whose
   endpoints are simple crossings. No parallels, no concurrency. This also
   kills the "at Tamura equality excess $=0$ forces full-side sharing" rescue:
   the line is tight yet double-covered.
2. **The extension of C–B Lemma 1 step (3) is false** ("each multiple point serves at most two
   shared-side pairs"): $n=6$, lines $y=0$, $y=x$, $y=-x$ (concurrent at the
   origin), $10y=10-x$, $x=-1$, $20y=x-20$; the six triangles
   $\{0,1,3\},\{1,2,3\},\{2,0,4\},\{0,1,4\},\{1,2,5\},\{2,0,5\}$ are pairwise
   interior-disjoint with **six** shared-side pairs at the single triple point.
3. **The per-line step fails under concurrency**: $n=6$, lines $y=0$,
   $y=x/5$, $x=1$, $10y=x-1$, $4y=x-2$, $3y=x-3$; the five pairwise-disjoint
   triangles $\{0,1,2\},\{0,3,4\},\{0,4,5\},\{0,1,3\},\{0,2,4\}$ all have a
   side on $y=0$, which carries only $4$ points / $3$ segments: $\sigma = 5 >
   n-2 = 4$. Thus the located segment-count argument does not establish
   $\lfloor n(n-2)/3\rfloor$ for the broad degenerate convention.

**What survives (proved by the same agent):**

* *Per-line lemma (no multipoint)*: putting line $\ell$ on the x-axis with
  crossing points $P_1 < \cdots < P_{k}$, the apex of the triangle on side
  $[P_i,P_j]$ ($g\ni P_i$, $h\ni P_j$) is above $\ell$ iff
  $\cot\theta_g > \cot\theta_h$; above-sides are strict cot-descents,
  below-sides strict ascents, and each same-side family is interior-disjoint,
  so $\sigma_\ell \le k_\ell - 1$. Consequently **if no three lines are
  concurrent, $3T \le n(n-2) - 2Q$**; at $n=10, T=26$: $Q \le 1$ — the
  has-parallel/no-concurrency cube has **exactly one parallel pair** and all
  other pairs crossing (`zparnoc`).
* *With multipoints*: $3T \le n(n-2) - 2Q - \sum_p k_p(k_p-4)$ (superseded by
  the full theorem in §7, which adds the crossing term $C$).
* *D₁₀ orbit transversal verified* (Burnside: $(120+5\cdot8)/20 = 8$ orbits;
  unit variables in the cube files decode to exactly the 8 representatives).
  Caveat: 8 is correct only while reflection remains unused; a second
  orientation-breaking clause would force the 12-orbit cyclic-only split.


## 7. The crossing-refined capacity theorem

Proved 2026-08-17 (independent construction by the proof agent, independent
adversarial audit `CrossingTheoremAudit`: verdict VALID, airtight; a second
exact adversarial audit `SuddenPtarmigan` attempted to break the token and
chord injections on 3,016 arrangements / 812,486 valid families /
5,900,863 linewise checks with planted pencils $k\in\{3,4,n{-}2,n{-}1\}$ and
shared-base collisions: **zero violations**, tight lines saturating
$\sigma+\gamma=m+2a$ 106,219 times; novelty
audit `scratch/kobon/novelty_audit_capacity_theorem.md`: no published source
contains the $C$ term or the $k_p(k_p-4)$ penalty in the broad convention).

**Theorem (capacity).** Let $A$ be an essential arrangement of $n$ distinct
affine lines ($Q$ = number of parallel pairs; finite points $p$ of
multiplicity $k_p\ge 3$). Let $S$ be any family of nondegenerate triangles
with sides on three lines of $A$ and pairwise disjoint **open** interiors —
crossed triangles, parallel pairs and multiple points allowed. Let $c_T$
count the lines of $A$ not supporting $T$ that meet the open interior of $T$,
and $C=\sum_T c_T$. Then

$$3|S| + C \;\le\; n(n-2) - 2Q - \sum_p k_p(k_p-4).$$

**Proof (gap-and-token injection).** Fix a line $\ell$ with $q_\ell$ parallel
mates, $\sigma_\ell$ selected sides on $\ell$, $\gamma_\ell$ selected
triangles whose open interior $\ell$ meets. In coordinates $\ell=\{v=0\}$,
write each line through a vertex $P_i=(u_i,0)$ as $u = u_i + c_g v$; set
$\alpha_i=\min c_g$, $\beta_i=\max c_g$. A base $[P_i,P_j]$ is above $\ell$
iff its endpoint slopes satisfy $c_g>c_h$, below iff $c_g<c_h$. Call an
elementary gap $e_t=[P_t,P_{t+1}]$ a *D-gap* if $\beta_t>\alpha_{t+1}$ and an
*A-gap* if $\alpha_t<\beta_{t+1}$. Chaining $\alpha\le\beta$ shows every
above-base contains a D-gap and every below-base an A-gap; map each base to
its leftmost such gap. Two above-bases cannot share a gap (their interiors
share a half-disk above $\ell$), nor two below-bases. An above/below collision
on $e_t$ forces a multipoint endpoint (else $\beta_t\le\alpha_{t+1}$ and
$\alpha_t\ge\beta_{t+1}$ contradict $c_t$ comparisons); keep the upper image
and send the lower image to a directional token of that multipoint (two
tokens per multipoint, reachable only from the adjacent gap). Selected sides
therefore inject into used gaps plus $2a_\ell$ tokens. If $\ell$ meets the
open interior of $T$, the open chord $T\cap\ell$ spans two distinct finite
vertices — including the $(0,+,-)$ corner case at a concurrent vertex — and
its leftmost interior elementary gap is unused by every selected side and
distinct for distinct $T$. Hence
$\sigma_\ell+\gamma_\ell \le m_\ell + 2a_\ell$ with $m_\ell$ = number of
bounded elementary gaps. Counting incidences,
$m_\ell+2a_\ell = n-2-q_\ell-\sum_{p\in\ell}(k_p-4)$; summing over $\ell$
gives the theorem (each $p$ lies on $k_p$ lines; $\sum q_\ell = 2Q$).
The all-parallel class is the exact exception ($|S|=C=0$ but RHS $=-n$).
$\blacksquare$

**Corollary (case trees).** At $(n,T)$ the slack is
$\Delta = n(n-2)-3T$:

* $n=12, T=39$: $\Delta=3$ — $C + 2Q + \sum k_p(k_p-4)\le 3$; without triple
  points $C\le 3$ ($Q{=}0$) resp. $C\le 1$ ($Q{=}1$).
* $n=14, T=54$ and $n=18, T=94$: $\Delta=6$ — same tree shape: $Q\le 3$, and
  $Q\ge 4$ needs triple points; each triple point adds exactly $3$ units of
  capacity; $k=4$ points are free; $k\ge 5$ costs $\ge 5$.
* $n=20, T=117$: $\Delta=9$ — $Q\le 4$.

The best known 14-line/53-triangle geometry (two independently verified
exact rational realizations with $Q\in\{2,3\}$ and no crossings:
`scratch/kobon/discoveries/n14_bader53_verified.json`,
`scratch/kobon/discoveries/n14_best_exact53.json`) sits at $3\cdot53=159$
against the $Q=3$ capacity ceiling $168-6=162=3\cdot54$: a 54-family inside
the no-concurrency $Q=3$ regime would saturate the capacity theorem with
equality on every line. Any 54-family must live at $Q\le 3$ with at most
simple concurrency patterns — exactly the cube tree under attack.

**Sharpness.** The bound is attained exactly by published record
arrangements at $n=9$ (Blanc-type $21$: $3\cdot21+0=63=9\cdot7$) and $n=15$
(rotationally symmetric $65$: $195=195$), and nearly tight on the certified
$n=10$ record ($75+0$ vs $78$, $Q=1$), the $n=11$ record ($96$ vs $99$) and
the $n=13$ record ($141$ vs $143$) — all computed by exact rational
`face_multiplicity_counts` on the gallery arrangements (full
maximum-independent-set selection, crossed triangles allowed). No
constant-term strengthening is possible in general.

**Theorem M (Q3 branch closed, $n=14$).** Any arrangement of 14 lines with
exactly three parallel pairs and no three concurrent carries at most 53
pairwise-interior-disjoint triangles. Capacity forces $C=0$ and equality
$\sigma_\ell+\gamma_\ell=n-2-q_\ell$ per line; the side/chord injections
become bijections, so every bounded elementary gap is a triangle side — a
cell — at equality. At each vertex the four-cell XOR forces exactly two
opposite triangular cells; pattern analysis kills mixed-extreme vertices and
forces $n$ double-extreme vertices, whose direction-adjacency wedge lemma
restricts them to cyclically adjacent slope classes. The resulting extremity
multigraph on the 11-class direction cycle must solve
$m_{i-1}+m_i=2s_i$ with $0\le m_i\le s_is_{i+1}$: **no integral solution
exists** for all 165 paired placements $s=(2,2,2,1^{8})$ (linear-system
certificate) and all 12 triad placements $s=(3,1^{11})$ of the even
12-class cycle (alternating-sum $\pm4$ inconsistency). Controls (Blanc-9,
rot-15, pentagon, hexagon at $n=6$!) stay feasible, so the obstruction is
genuinely $n{=}14$-specific. Standalone 0.4 s stdlib certificate:
`scratch/kobon/n14/escape/appendix_independent_check.py`; snapped Bader
realization has exactly 3 escapes, so the branch value equals 53.
**Consequence:** the no-concurrency $Q=3$ part of the $K(14)$ case tree is
closed analytically — the live branches reduce to $Q\in\{0,1,2\}$ plus
concurrency signatures. The same even-cycle mechanism closes the $Q=3$
branches of $n{=}18$, $T{=}94$.

**Encoding.** `add_face_bound(..., crossing_lits=<exact TC family>,
per_line=True)` adds the per-line inequalities
$\sigma_r+\gamma_r \le n-2-q_r-\sum_{p\in r}(k_p-4)$; the no-concurrency
variant with constant $q_r$ is `add_no_concurrency_crossing_budget`. Both are
pure necessary-condition accelerators (§3.4 soundness architecture applies).

## 8. Artifacts

Canonical (this directory, `math/kobon/`):

- `engine.py` — full validated engine + reproducible instance driver
  (`dump_instances(10, 26)` regenerates all 12 proof-family CNFs
  byte-identically; sealed by the 2026-08-15 audit in a fresh venv).
- `report.md`, `AUDIT.md`, `README.md` — this report, the certification
  audit, and the SSOT status index.

Working data (`~/jinleic-workspace/scratch/kobon/`):

- `kobon_n10_t26_proof{,_simple,_par_noc,_par_conc,_c0..c7}.cnf` — the
  proof-relevant family (base + 11 exhaustive cubes; single WLOG σ clause,
  justified by $b \mapsto -b$).
- `proof_n10_t26_*.drat` — 9 complete cube proofs (empty clause at EOF) and
  the 2 truncated ones (`c0`, `par_conc`) being regenerated.
- `kobon_n{9,10,12,14}_t*.cnf`, `kobon_n10_t26_{exact2,x*,z*,s*,b*,…}.cnf` —
  historical instances from earlier engine revisions (sound but superseded).
- `proof_n9_t22.drat` — certified control: `s VERIFIED` by `drat-trim`.

Audit working dir (`~/jinleic-workspace/scratch/kobon-audit/`): tools
(kissat 4.0.4, drat-trim), `split/` subcube CNFs, `resolve/` replacement
proofs, `results.tsv` fleet verdicts, `memwatch.log`.
