# The Q=3, no-concurrency branch of K(18), target T=94: Theorem-M obstruction

**Verdict (conditional only on the capacity theorem and the geometric lemmas below):
OBSTRUCTION.** No arrangement of 18 distinct affine lines with exactly three
parallel pairs (in the paired layout) or one parallel triad (in the triad
layout), with no three concurrent, can carry 94 pairwise-open-interior-disjoint
triangles.

This file is an analytic derivation: the capacity-theorem hypothesis and the
geometric lemmas (E0–L4 and the proof-interface steps) are proved by hand in
§3–§6 and are **not** machine-checked. The finite cycle-enumeration legs
(the 455 paired placements, the 16 triad placements, the exact 19-orbit
multiplicity table, the 231 and 9 reflection quotients, and the controls)
ARE machine-checked by `scratch/kobon/n18/appendix_check_n18.py` (solver
functions copied verbatim from the audited n=14 appendix), executed
2026-08-20: **finite enumeration items VERIFIED**.

---

## 1. Setup and the equality regime

Use the capacity theorem in the form

$$
3|S|+C
\;\le\;
 n(n-2)-2Q-\sum_p k_p(k_p-4),
\qquad
 C=\sum_{T\in S}c_T,
$$

where $c_T$ is the number of non-supporting lines meeting the open interior of
$T$. In the no-concurrency branch there are no $k_p\ge 3$ points. For
$n=18$, $Q=3$, and $|S|=94$,

$$
3|S|+C=3\cdot94+C=282+C,
$$
$$
18(18-2)-2\cdot3=18\cdot16-6=288-6=282,
$$
so
$$
282+C\le282 \quad\Longrightarrow\quad C=0.
$$

The per-line capacity is $n-2-q_\ell$, where $q_\ell$ is the number of lines
parallel to $\ell$.

* **Paired layout:** six lines lie in three size-two parallel classes. Each
  paired line has $q_\ell=1$ and hence $15$ bounded gaps; the other twelve
  singleton lines have $q_\ell=0$ and hence $16$ bounded gaps. Thus
  $$6\cdot15+12\cdot16=90+192=282.$$
* **Triad layout:** three lines lie in one size-three parallel class. Each of
  those three lines has $q_\ell=2$ and hence $14$ bounded gaps; the other
  fifteen singleton lines have $16$ bounded gaps. Thus
  $$3\cdot14+15\cdot16=42+240=282.$$

Since $C=0$ and $3|S|=282$, equality holds in every per-line inequality. In
particular, every non-supporting-interior count is zero and every bounded
elementary gap is used by exactly one selected side. We call this equality
hypothesis

> **(SAT)** Every bounded elementary gap of the arrangement is the image of
> exactly one selected side.

The rest of the proof assumes (SAT) and derives a contradiction in both
parallel-class layouts.

### Lemma E0 (cells)

Under (SAT), every selected side is one elementary gap and every selected
triangle is an elementary cell (an empty zigzag triangle).

*Proof.* If a selected side on $\ell$ ran through an interior crossing
$P_k$ of $\ell$, the line through $P_k$ would enter the open interior of the
triangle, contradicting $C=0$. Hence all three selected sides are gaps. A
triangle bounded by gap segments and containing no line in its open interior is
exactly a face of the arrangement. $\square$

### Lemma E1 (one triangular user per gap)

Under (SAT), each gap has exactly one adjacent triangular cell. The other
adjacent cell is not a triangle.

*Proof.* Put $\ell$ horizontally and let $c(P_i)$ be the cotangent of the
crossing line at the $i$th crossing of $\ell$. For the gap
$e_t=[P_tP_{t+1}]$, the cell above it is triangular exactly when
$c_t>c_{t+1}$ and its apex triangle is uncrossed; the corresponding condition
below is $c_t<c_{t+1}$ and uncrossed. In the simple (no-concurrency) regime
$c_t\ne c_{t+1}$, so the two strict cot comparisons are mutually exclusive
(a flat gap is treated in Lemma L5). By (SAT) the gap is used once, and by E0
its user is a triangular cell. Therefore exactly one of the two adjacent cells
is that user. $\square$

There are $282$ used gaps and each selected cell uses three gaps, so the
number of selected cells is $282/3=94$. In particular, there is no escape in
this equality regime.

---

## 2. Vertex consequences of saturation

### Lemma L2 (zigzag alternation)

At every interior vertex of every line, the cotangent sequence alternates
strictly:

$$
c_{t-1}<c_t>c_{t+1}
\quad\text{or}\quad
c_{t-1}>c_t<c_{t+1}.
$$

*Proof.* At $P=\ell\cap f$, the four local gaps ($e^-,e^+$ on $\ell$ and
$g^-,g^+$ on $f$) bound four quadrant cells $x_1,x_2,x_3,x_4$. By E1 each
bounded gap has exactly one triangular adjacent cell, and a triangular
quadrant cell covers exactly the two gaps bounding that quadrant. If
$t_j\in\{0,1\}$ records whether $x_j$ is triangular, the four gap equations
are

$$
t_1+t_2=t_2+t_3=t_3+t_4=t_4+t_1=1.
$$

Their only solutions are $(1,0,1,0)$ and $(0,1,0,1)$: two opposite quadrant
cells are triangles. Reading the two possible choices with the cot/side
convention gives an ascent followed by a descent or a descent followed by an
ascent. The inequalities are strict because there is no concurrency. $\square$

### Lemma L3 (no mixed extreme; exactly 18 double-extreme vertices)

Call $P=\ell\cap f$ extreme on $\ell$ if it is the first or last crossing on
$\ell$, and double-extreme if it is extreme on both lines. Under (SAT):

1. no vertex is extreme on exactly one of its two lines; and
2. the number of double-extreme vertices is $18$.

*Proof.* Suppose $P$ is rightmost on $f$ but interior on $\ell$. The two cells
adjacent to the empty ray of $f$ at $P$ are unbounded, hence cannot be
triangles. Both gaps of $\ell$ incident to $P$ are bounded. By E1 their
triangular users must be the other cells adjacent to those gaps; both of
those cells are also adjacent to the bounded gap of $f$, contradicting E1's
one-sidedness. The other orientations are identical, proving (1).

There are exactly two extreme incidences on each of the 18 lines, hence
$2\cdot18=36$ extreme incidences. If $c$ is the number of double-extreme
vertices and $b$ the number of mixed-extreme vertices, then
$36=2c+b$. Part (1) gives $b=0$, so $c=18$. $\square$

### Lemma L4 (direction adjacency)

At every double-extreme vertex, the two support directions are adjacent in
the circular order of slope classes: no third direction lies strictly inside
the cone bounded by the two empty rays.

*Proof.* The wedge $W$ between the two empty rays is a convex open cone of
angle less than $\pi$. A line $g$ whose direction is strictly inside this cone
crosses both support lines. Since the crossing is not on an empty ray, both
crossings would have to lie on the non-empty rays. In cone coordinates, a line
through points $(-\alpha,0)$ and $(0,-\beta)$ with $\alpha,\beta>0$ misses the
open positive quadrant, so it cannot enter $W$. This contradicts that its
direction lies strictly inside the wedge. Therefore the two support directions
are adjacent. $\square$

### Lemma L5 (flat gaps cannot occur under SAT)

A flat gap whose endpoints are crossings with the two members of one parallel
pair and which contains no other crossing has two unbounded half-strip cells
adjacent to it, so it cannot be used. Hence (SAT) rules out every flat gap.

*Proof.* The two parallel lines do not meet, and with no third crossing between
their intersections on $\ell$, the two cells on the sides of the gap continue
along the corresponding empty rays. Both are unbounded. E1 says a used gap has
a triangular user, a contradiction. $\square$

### Lemma L6 (a gap is not both a descent and an ascent)

In the simple regime, no elementary gap is simultaneously a descent and an
ascent. At either endpoint only one line crosses $\ell$, so the cot comparison
is total: a gap is a strict descent, a strict ascent, or a flat gap. By L5 the
flat case is absent under (SAT).

---

## 3. The extremity multigraph and the odd $R=15$ leg

### Theorem M18

No simple (no-three-concurrent) arrangement of 18 lines with $Q=3$ satisfies
(SAT). Consequently, the no-concurrency $Q=3$ branch cannot attain $T=94$.

*Proof.* By L3 there are 18 double-extreme vertices. By L4 each is at the
crossing of two cyclically adjacent slope classes. Let
$C_1,\ldots,C_R$ be the slope classes in circular order, with sizes $s_i$.
Define

$$
m_i=\#\{\text{double-extreme vertices at }C_i\cap C_{i+1}\},
\qquad i\in\mathbb Z_R.
$$

Every line in $C_i$ contributes its two extreme vertices, both double-extreme
and each between adjacent classes. Therefore

$$
\tag{M1}m_{i-1}+m_i=2s_i.
$$

There are only $s_is_{i+1}$ crossings between two neighboring classes, and
different double-extreme vertices use different line pairs. Therefore

$$
\tag{M2}0\le m_i\le s_is_{i+1},
\qquad m_i\in\mathbb Z.
$$

We now treat the two possible $Q=3$ class structures.

### 3A. Paired layout: $(2,2,2)+12$ singles, $R=15$ (odd)

The three size-two classes give three doubled positions among the 15 positions
of the direction cycle. Let their cyclic positions be $p_1,p_2,p_3$, and let
$r_j\ge0$ be the number of singleton classes strictly between $p_j$ and
$p_{j+1}$, cyclically. Then

$$r_1+r_2+r_3=15-3=12.$$

Set

$$x_i=s_i-1\in\{0,1\},\qquad y_i=m_i-1.$$

Then (M1) becomes

$$y_{i-1}+y_i=2x_i.\tag{3.1}$$

Let $a_j=y_{p_j}$ be the value on the edge leaving the $j$th doubled
position. Across the $r_j$ singleton equations, the sign alternates, and the
next doubled equation gives the exact recurrence

$$a_{j+1}=2-(-1)^{r_j}a_j.\tag{3.2}$$

Because $r_1+r_2+r_3=12$ is even, exactly one of the following parity cases
occurs.

* **All three gaps even.** Equation (3.2) has the unique solution
  $a_1=a_2=a_3=1$. Along every gap, $y$ is therefore
  $$1,-1,1,-1,\ldots,$$
  so any gap of length at least $4$ contains a single-single edge with
  $y=1$, i.e. $m=2$. On a single-single edge (M2) has upper bound
  $$s_is_{i+1}=1,$$
  so this is impossible. If all three even gaps had length at most $2$, their
  sum would be at most $6$, contrary to $r_1+r_2+r_3=12$. Thus every all-even
  placement is infeasible by a bound violation.

* **Exactly two gaps odd and one gap even.** Relabel cyclically so the parity
  sequence is odd, odd, even. Solving (3.2) gives, up to cyclic relabeling,
  $$
  (a_1,a_2,a_3)=(-1,1,3),
  $$
  hence the three edges leaving doubled classes have
  $$m=(0,2,4).
  $$
  The $m=4$ edge is the first edge of the even gap. If that even gap has
  length $e>0$, this edge is doubled-single and has upper bound $2$, so
  $m=4>2$; the next edge has $y=-3$, hence $m=-2<0$. If $e=0$, the $m=4$
  edge is doubled-doubled and its upper bound $4$ is not itself a violation.
  However, an odd gap of length at least $3$ contains an alternating
  single-single edge with $y=1$, hence $m=2>1$. Avoiding that violation in
  both odd gaps would force both to have length $1$, which would make the
  total gap sum $1+1+0=2$, not $12$. Therefore the $e=0$ placements are
  also infeasible.

The two cases exhaust all placements, so the paired layout has no
(M1)+(M2)-feasible solution. The exact placement counts and the reflection
arithmetic are recorded next.

#### Exact placement and orbit count

There are

$$\binom{15}{3}=455$$

raw placements of three indistinguishable doubled classes among 15 positions.
Equivalently, a distinguished doubled position and an ordered gap composition
can be counted as

$$
\frac{15\binom{12+3-1}{3-1}}{3}
=\frac{15\binom{14}{2}}3
=\frac{15\cdot91}{3}=455,
$$

because every unrooted placement has three choices of distinguished doubled
position.

For a fixed reflection of a 15-cycle, there is one fixed position and seven
transposed pairs. A reflection-invariant 3-subset must contain the fixed
position and one transposed pair, giving $7$ fixed placements. Thus the
reflection-only quotient has exactly

$$\frac{455+7}{2}=231$$

orbits.

If positions are regarded as circularly relabeled as well (the useful orbit
notion for the gap table), Burnside over the dihedral group $D_{15}$ gives
exactly 19 orbits. The identity rotation fixes $455$ placements; the two
rotations by $5$ and $10$ positions each fix $5$ placements (one 3-cycle of
positions); the other twelve nonidentity rotations fix none. Each of the 15
reflections fixes $7$ placements. Hence

$$
\frac{455+2\cdot5+15\cdot7}{30}
=\frac{570}{30}=19.
$$

These 19 circular/reflection orbit representatives are the unordered gap
triples $\{r_1,r_2,r_3\}$, i.e. the partitions of 12 into at most three
nonnegative parts. A representative with all three entries distinct accounts
for $5\cdot6=30$ raw placements; one with exactly two equal entries accounts
for $5\cdot3=15$; and $(4,4,4)$ accounts for $5$.

| gap representative $\{r_1,r_2,r_3\}$ | parity / even-gap case | raw placements | exact obstruction |
|---|---:|---:|---|
| $(0,0,12)$ | all even | 15 | $m=2$ on a single-single edge |
| $(0,1,11)$ | two odd, $e=0$ | 30 | $m=2>1$ on a single-single edge |
| $(0,2,10)$ | all even | 30 | $m=2$ on a single-single edge |
| $(0,3,9)$ | two odd, $e=0$ | 30 | $m=2>1$ on a single-single edge |
| $(0,4,8)$ | all even | 30 | $m=2$ on a single-single edge |
| $(0,5,7)$ | two odd, $e=0$ | 30 | $m=2>1$ on a single-single edge |
| $(0,6,6)$ | all even | 15 | $m=2$ on a single-single edge |
| $(1,1,10)$ | two odd, $e=10>0$ | 15 | $m=4>2$ and $m=-2<0$ |
| $(1,2,9)$ | two odd, $e=2>0$ | 30 | $m=4>2$ and $m=-2<0$ |
| $(1,3,8)$ | two odd, $e=8>0$ | 30 | $m=4>2$ and $m=-2<0$ |
| $(1,4,7)$ | two odd, $e=4>0$ | 30 | $m=4>2$ and $m=-2<0$ |
| $(1,5,6)$ | two odd, $e=6>0$ | 30 | $m=4>2$ and $m=-2<0$ |
| $(2,2,8)$ | all even | 15 | $m=2$ on a single-single edge |
| $(2,3,7)$ | two odd, $e=2>0$ | 30 | $m=4>2$ and $m=-2<0$ |
| $(2,4,6)$ | all even | 30 | $m=2$ on a single-single edge |
| $(2,5,5)$ | two odd, $e=2>0$ | 15 | $m=4>2$ and $m=-2<0$ |
| $(3,3,6)$ | two odd, $e=6>0$ | 15 | $m=4>2$ and $m=-2<0$ |
| $(3,4,5)$ | two odd, $e=4>0$ | 30 | $m=4>2$ and $m=-2<0$ |
| $(4,4,4)$ | all even | 5 | $m=2$ on a single-single edge |
| **total** | **7 all-even; 3 with $e=0$; 9 with $e>0$** | **455** | **all infeasible** |

The raw-count subtotal is exact:

$$
\underbrace{140}_{\text{all-even}}
+\underbrace{90}_{\text{two odd},\ e=0}
+\underbrace{225}_{\text{two odd},\ e>0}
=455.
$$

Thus the odd $R=15$ system has a unique candidate for every raw placement,
and all $455$ candidates violate (M2) (by a negative entry and/or an upper
bound). The 19 gap representatives above enumerate those placements modulo
circular relabeling and reflection; the reflection-only count is 231 as shown
above.

### 3B. Triad layout: one size-three class + 15 singles, $R=16$ (even)

For the actual 18-line triad layout, one slope class has size $3$ and the
remaining $15$ classes are singletons:

$$3+15=18,\qquad R=1+15=16.
$$

There are 16 raw placements of the size-three class around the even direction
cycle. The necessary alternating-sum condition follows directly from (M1):

$$
0=\sum_{i=0}^{15}(-1)^i(m_{i-1}+m_i)
 =\sum_{i=0}^{15}(-1)^i(2s_i).
$$

For a placement at position $t$, the all-single baseline cancels because
$\sum_{i=0}^{15}(-1)^i=0$. Replacing $s_t=1$ by $s_t=3$ adds $4$ at position
$t$, so

$$
\sum_{i=0}^{15}(-1)^i(2s_i)=4(-1)^t\in\{+4,-4\}\ne0.
$$

Therefore no placement has even a rational solution of (M1), let alone an
integer solution satisfying (M2): exactly $16/16$ triad placements are
infeasible. For a fixed reflection of the 16-cycle there are two fixed
positions and seven transposed pairs, so the reflection-only placement count is
$(16+2)/2=9$ orbits; rotations make all 16 positions one dihedral orbit.

This proves Theorem M18. $\square$

---

## 4. Consequence for the no-concurrency branch

A 94-family at $(n,Q)=(18,3)$ would force (SAT) by the capacity equality.
Theorem M18 rules out (SAT) in both the three-pair and one-triad layouts.
Hence the no-concurrency $Q=3$ branch cannot attain 94; its capacity-theorem
upper bound is sharpened to

$$|S|\le93.$$

This conclusion uses no numerical arrangement search. It is an analytic
obstruction once the generic capacity theorem and the geometric saturation
lemmas are accepted.

---

## 5. Exact arithmetic ledger

For later checking, the complete hand arithmetic is:

* $18\cdot16=288$;
* $2Q=2\cdot3=6$;
* $288-6=282$;
* $3T=3\cdot94=282$;
* paired gaps: $6\cdot(18-2-1)+12\cdot(18-2)=6\cdot15+12\cdot16=282$;
* triad gaps: $3\cdot(18-2-2)+15\cdot(18-2)=3\cdot14+15\cdot16=282$;
* paired slope classes: $3\cdot2+12\cdot1=18$, $R=3+12=15$;
* triad slope classes: $1\cdot3+15\cdot1=18$, $R=1+15=16$;
* paired placement count: $\binom{15}{3}=455$;
* reflection-only paired quotient: $(455+7)/2=231$;
* dihedral paired quotient: $(455+2\cdot5+15\cdot7)/30=19$;
* paired obstruction subtotals: $140+90+225=455$;
* triad placements: $16$, alternating sum $4(-1)^t=\pm4$.

---

## 6. Hypothesis ledger

| item | hypothesis | conclusion | status in this hand artifact |
|---|---|---|---|
| Capacity equality | capacity theorem, $n=18,Q=3,T=94$, no multipoints | $C=0$, all linewise bounds tight | algebra shown |
| E0/E1 | (SAT) and simple arrangement | sides are gaps; one triangular user per gap | proof copied/adapted from the gap geometry |
| L2 | (SAT), interior simple vertex | strict zigzag alternation | local XOR proof |
| L3 | (SAT) | no mixed extreme; $c=18$ | extremity count proof |
| L4 | simple arrangement, double-extreme vertex | adjacent slope directions | wedge proof |
| L5/L6 | simple arrangement (and SAT for L5) | no flat or both-D/A gap | local gap proof |
| M paired leg | (M1)+(M2), $R=15$ | all 455 placements infeasible | hand algebra; exact check pending |
| M triad leg | (M1), $R=16$ | all 16 placements have alternating sum $\pm4$ | hand algebra; exact check pending |

---

## 7. Verification queue — **ALL ITEMS VERIFIED 2026-08-20**

The queue was executed by `scratch/kobon/n18/appendix_check_n18.py` (pure
stdlib, exact integers; solver functions copied verbatim from the audited
`scratch/kobon/n14/escape/appendix_independent_check.py`). Observed output
(`python3 scratch/kobon/n18/appendix_check_n18.py`):

1. **VERIFIED ODD-R15 / 455 placements.** All $\binom{15}{3}=455$ paired
   placements visited; every one ruled infeasible by the extracted
   bound-check (`455/455 paired placements infeasible`). Hand subtotals
   reproduced exactly: `{'all_even': 140, 'two_odd_e0': 90,
   'two_odd_ep': 225}`. `19 dihedral orbit representatives confirmed`;
   reflection-fixed placements = 7 and the reflection-only quotient
   $(455+7)/2=231$ confirmed.
2. **VERIFIED EVEN-R16 / triad.** All 16 placements infeasible; the
   alternating sum equals $4(-1)^t$ (checker pattern `+-+-+-+-...`,
   $|a|=4$ everywhere); the reflection-only quotient $(16+2)/2=9$
   confirmed.
3. **VERIFIED controls and proof-interface checks.** All-single odd/even
   cycles (`n=9`, `R=15`, `R=17`, pentagon, even `R=4`) and the $(2,2,2)$
   hexagon control all return feasible, so the checker is not
   over-strong.
