# A Crossing-Refined Capacity Theorem for Families of Triangle Selections in Arrangements of Lines

**Companion to** `math/kobon/report.md` (the certified $K_{\rm gen}(10)=25$ decision). This
document is the standalone, publication-grade writeup of the capacity theorem announced in
report.md §7, including the full proof, its sharpness analysis, the comparison with the
prior literature, the consequences used to drive the $n\in\{11,12,14,18,20\}$
decision trees, and the computational validation behind it.

What is *not* claimed here: no new exact value of the Kobon number is decided by this
theorem alone; in particular **$K(11)$, $K(12)$, $K(14)$ are not settled here**
(the companion SAT/DRAT campaigns are open, see §8). The theorem is a *necessary
condition*: it bounds every admissible selection in the broad convention, and thereby
carves the target-$T$ search space in a way no published bound does.

---

## Abstract

Let $\mathcal{A}$ be an essential arrangement of $n$ distinct affine lines (at
least one nonparallel pair), and let $S$ be any family of nondegenerate triangles
with sides on $\mathcal{A}$ and pairwise disjoint **open** interiors — other lines may
cross a selected triangle, parallel lines and multiple points are allowed. With $Q$ the
number of parallel pairs, $k_p\ge3$ the multiplicities of the finite multiple points,
and $C=\sum_{T}\#\{\ell\not\ni T:\ \ell\cap\operatorname{int}(T)\neq\emptyset\}$
the total number of (triangle, outside line) interior crossings, we prove

> **Theorem.** $3|S|+C\;\le\; n(n-2)-2Q-\sum_p k_p(k_p-4).\tag{1}$

The proof is elementary and needs no general-position hypothesis. Fix a line $\ell$. A
"gap-and-token" injection maps the selected sides on $\ell$ injectively into the bounded
elementary gaps of $\ell$ plus two formal directional tokens at each multiple point (an
above/below collision forces a multiple endpoint, which absorbs it); a "chord" injection
maps each triangle whose open interior $\ell$ meets to a further vacant elementary gap.
Thus $\sigma_\ell+\gamma_\ell\le m_\ell+2a_\ell
=n-2-q_\ell-\sum_{p\,\text{on}\,\ell}(k_p-4)$, and summing over the lines yields
(1). The all-parallel class is the unique exception (essentiality is the exact
hypothesis; there the RHS is $-n$).

The bound is sharp as stated: equality holds at $n=9$ and $n=15$; and tight within
$3,3,2$ units at $n=10,11,13$ — each value computed exactly
(`face_multiplicity_counts`; independently reproduced for this paper). A novelty audit
(`scratch/kobon/novelty_audit_capacity_theorem.md`) found **no published source**
containing $C$ or $\sum_p k_p(k_p-4)$, and exact-rational witnesses show the local
charging steps of the classical arguments (Tamura; the Clément–Bader draft; BBL
arXiv:0706.0723; Blanc arXiv:0801.2845) fail precisely in this convention. An
independent proof audit (`CrossingTheoremAudit`) returned **VALID** and re-derived (1);
a further adversarial audit (`SuddenPtarmigan`) executed **$5{,}900{,}863$ linewise
exact checks** over $3{,}016$ arrangements and $812{,}486$ valid selections with
**zero violations**. The inequality is then used to derive the audited budget case trees
that drive the $n\in\{11,12,14,18,20\}$ decision campaign; those values remain open.

---

## 1. Conventions

**Arrangements.** An *arrangement of lines* is a finite set of distinct affine lines
$\ell:y=mx+b$ in $\mathbb{R}^2$. All arrangements here are *essential*: not all
lines are parallel. Write

- $Q$ for the number of parallel pairs $\{\ell,\ell'\}$;
- $q_\ell$ for the number of lines parallel to $\ell$.

For a line $\ell$, the *finite vertices* on $\ell$ are the points of intersection
of $\ell$ with some other line. Because the arrangement is essential, every line carries
at least one finite vertex (§4, Lemma 1). A finite point $p$ at which exactly
$k_p$ lines meet is a *multiple point* of multiplicity $k_p\ge 2$; write
$k_p\ge3$ when the point is a genuine *multipoint*, and let $a_\ell$ be the number
of multipoints on $\ell$. Between consecutive finite vertices $P_t,P_{t+1}$ of
$\ell$ ($P_1<\cdots<P_r$ in the natural order along $\ell$) lies the *elementary
gap* $e_t=[P_t,P_{t+1}]$; the number of bounded elementary gaps on $\ell$ is
$m_\ell=r-1$.

**Selections.** A *triangle* is nondegenerate if its three vertices are not collinear.
A *selected* triangle has each of its three sides on a line of $\mathcal{A}$ (two
sides may the same line only in a degenerate configuration, which we exclude; the three
sides lie on three distinct lines, and a triangle is *crossed* by any arrangement line
meeting its open interior). $S$ is a family of nondegenerate triangles with pairwise
disjoint **open** interiors; the closures may overlap, share boundary segments, or meet at
points and segments. This is the *broad convention* of report.md §1, in which
$K_{\rm gen}(n)$ denotes the largest possible $|S|$. For $T\in S$ and a line
$\ell\in\mathcal{A}$ not supporting $T$, define $c_T(\ell)=1$ if
$\ell\cap\operatorname{int}(T)\neq\emptyset$ and $0$ otherwise; let $c_T$ be the
number of $\ell$ with $c_T(\ell)=1$, and $C=\sum_{T\in S}c_T$. Finally let $\sigma_\ell=\#\{T\in S: \ell$ supports a side
of $T\}$ and $\gamma_\ell=\#\{T\in S: \ell\cap\operatorname{int}(T)\neq\emptyset\}$.
Then $\sum_\ell \sigma_\ell=3|S|$ and $\sum_\ell\gamma_\ell=C$.

> **Remark (why $C$ matters).** In the classical *face* convention one insists that a
> selected triangle is a triangular face of the arrangement; then $C\equiv 0$ by
> definition and $|S|=T\le$ (number of bounded triangular faces). The broad
> convention is strictly more general: report.md §2 records
> $K_{\rm gen}(8)\ge 15 > as_3(8)=14$ and $K_{\rm gen}(12)\ge 38>as_3(12)=37$,
> so crossed triangles genuinely beat the face count. Condition (1) is the first bound
> that prices them.

---

## 2. The theorem

**Theorem 1** (capacity, crossing-refined). *Let $\mathcal{A}$ be an essential
arrangement of $n$ distinct affine lines, let $S$ be any family of nondegenerate
triangles with sides on $\mathcal{A}$ and pairwise disjoint open interiors, and let
$C$ be as above. Then*

$$3|S|+C\;\le\; n(n-2)-2Q-\sum_p k_p(k_p-4),\tag{1}$$

*where $p$ ranges over the finite multiple points of $\mathcal{A}$ with $k_p\ge 3$.*

**Corollary 2** (integer slack form). At a target $T$ the *slack* is
$\Delta=n(n-2)-3T$ and every $T$-family must satisfy

$$C+2Q+\sum_p k_p(k_p-4)\le \Delta.\tag{2}$$

In particular, if $T$-families exist, $2Q+\sum_p k_p(k_p-4)<\Delta$ rules them
out when $C$ forces a lower bound, and any exterior crossing of a selected triangle
must be "paid for" out of $\Delta$.

**Corollary 3** (simple case). If $\mathcal{A}$ is simple (no parallel pair, no
multipoint), then $3T\le n(n-2)$, i.e. $T\le\lfloor n(n-2)/3\rfloor$. This
recovers Tamura's bound in the simple convention; the novelty of (1) is the
degeneracy- and crossing-refined right-hand side. [Simple special case is
folklore/easily derivable; see novelty audit §Verdict(4b).]

The statement of Theorem 1 matches exactly the statement proven in `report.md` §7;
the proof below is the gap-and-token + chord proof from §7, expanded to full
formality following the audited version.

---

## 3. The proof

The idea. Fix a line $\ell$. We exhibit two disjoint injections whose codomains
together have size $n-2-q_\ell+2a_\ell$: (i) an injection of the selected
*side-intervals* on $\ell$ (there are $\sigma_\ell$ of them) into a set of
"occupied positions" consisting of the $m_\ell$ elementary gaps plus $2a_\ell$
directional tokens at multipoints, and (ii) an injection of the selected *crossings*
(counted by $\gamma_\ell$) into the elementary gaps not used by (i). Then
$\sigma_\ell+\gamma_\ell\le m_\ell+2a_\ell$, and summing gives (1).

No general-position assumption is needed anywhere: the whole point of the token set is to
absorb the collisions that multiple points create, and the chord map is forced to use
*unused* gaps, so it never collides with the side map.

### 3.1 Setup on one line

Fix $\ell$, and put it in coordinates $\ell=\{v=0\}$. Let
$P_1<\cdots<P_r$ be its finite vertices, ordered along $\ell$. Every line $g\ne\ell$
meeting $\ell$ at $P_t$ is, in these coordinates, a graph $u=u_t+c_g v$ (this is
the $x=x_i+c_i y$ normalization of report.md §3.5; $c_g$ is the reciprocal
slope convention). For $P_t$ write

$$\alpha_t=\min\{c_g:\ g\ni P_t,\ g\ne\ell\},\qquad \beta_t=\max\{c_g:\ g\ni P_t\},\qquad \alpha_t\le\beta_t,$$

with $\alpha_t=\beta_t$ exactly when $P_t$ is simple.

**Base criterion.** Let $T$ have the side $[P_i,P_j]$ ($i<j$) on $\ell$, with the
apex $A= g\cap h$, $g\ni P_i$, $h\ni P_j$. In these coordinates $T$ lies
*above* $\ell$ (i.e. $v>0$ near $[P_i,P_j]$ in $T$) exactly when $c_g>c_h$,
and *below* exactly when $c_g<c_h$. This is the elementary lemma behind report.md
§3.5's "no-concurrency side capacity": it is a strict comparison of the two slopes and
holds verbatim at multiple vertices.

### 3.2 The side map: gaps and tokens

Call the elementary gap $e_t$ a **D-gap** if $\beta_t>\alpha_{t+1}$ and an
**A-gap** if $\alpha_t<\beta_{t+1}$.

**Lemma 1.** *An above-base contains an elementary D-gap in its relative interior;
a below-base contains an A-gap in its relative interior.*

*Proof.* Let $[P_i,P_j]$, $i<j$, be an above-base, $c_g>c_h$. If no $t$,
$i\le t<j$, is a D-gap, then $\beta_t\le \alpha_{t+1}$ for every such $t$, and
chaining with $\alpha_s\le\beta_s$ gives
$c_g\le \alpha_{i+1}\le\beta_{i+1}\le\cdots\le\alpha_j\le c_h$, contradicting
$c_g>c_h$. The below-case chains the reverse inequalities. Since $[P_i,P_j]$
is the intersection of the half-planes of its endpoints, the found gap lies in its
relative interior. ∎

**Lemma 2** (same-side injectivity). *Two above-bases on $\ell$ cannot contain the
same elementary gap; neither can two below-bases.* Equivalently, two same-side
selected sides have relative-interior-disjoint base intervals (they may meet only at
endpoints).

*Proof.* All base endpoints are arrangement vertices, so a positive-length overlap of two
base intervals contains the relative interior of an elementary gap $e$. At a point
$x\in\operatorname{relint}(e)$ each of the two triangles contains a small open half-disk
above (resp. below) $\ell$ in its interior for sufficiently small radius, so the two
open interiors intersect, contradicting disjointness. Bases meeting only at an endpoint
share no elementary gap. ∎

Map each above-base to its **leftmost** D-gap and each below-base to its leftmost
A-gap. By Lemma 2 the maps restricted to above-bases resp. below-bases are injective.

**Lemma 3** (collisions need multipoints; token resolution). *If an above-base and a
below-base are both mapped to the same elementary gap $e_t$, then $e_t$ is both a
D-gap and an A-gap, and at least one of $P_t,P_{t+1}$ is a multipoint.*

*Proof.* The gap is simultaneously D and A by the two images. If both endpoints were
simple then $\alpha_t=\beta_t=c_t$, $\alpha_{t+1}=\beta_{t+1}=c_{t+1}$, and
"D" says $c_t>c_{t+1}$ while "A" says $c_t<c_{t+1}$: contradiction. ∎

**Important scope note** (audited). The conclusion of Lemma 3 applies to an *image
collision*, not to physical overlapping segments. A shared elementary segment need **not**
have a multiple endpoint — this is exactly failure mode (1) of report.md §6, and no
strengthening to "shared segment $\Rightarrow$ multipoint" is claimed (the exact
counterexample is in §5.3 below).

Now fix the image collision. For each multipoint $P_t$ on $\ell$ introduce two formal
*tokens* $P_t^-$ (reachable only from the left adjacent gap $e_{t-1}$) and $P_t^+$
(reachable only from $e_t$). When a (unique, by Lemma 2) lower image collides with an
upper image on $e_t$, retain the upper image on the gap: send the lower image to $P_t^+$
if $P_t$ is a multipoint, otherwise to $P_{t+1}^-$.

**Lemma 4** (token injectivity). *The full side map — gaps for first images, tokens
for displaced lower images — is injective and uses only $2a_\ell$ tokens.* In
particular every selected side on $\ell$ is accounted for, and
$\sigma_\ell\le m_\ell+2a_\ell$.

*Proof.* A collision gap $e_t$ carries exactly one lower image (two lower images would
share a gap, forbidden by Lemma 2), so at most one displacement occurs per gap. The
token $P_t^+$ is reachable only from $e_t$, $P_t^-$ only from $e_{t-1}$, and
adjacent collision gaps $e_{t-1},e_t$ therefore use $P_t^-,P_t^+$
— two distinct formal positions; no token is reused. Endpoints at the ends of the vertex
chain have at most one usable direction; an outward token is simply unused. Thus the image
set lives in $G_\ell\cup\{P_t^-,P_t^+:\ P_t\ \text{a multipoint}\}$, of size
$m_\ell+2a_\ell$. The apex-over-apex case (two above-bases) can never collide by
Lemma 2. ∎

### 3.3 The chord map: crossings occupy unused gaps

Let $\ell$ cross the open interior of a selected triangle $T$ (so $\ell$ does not
support $T$).

**Lemma 5.** *$\ell\cap\operatorname{int}(T)$ is a nonempty open chord whose two
endpoints are distinct finite arrangement vertices.*

*Proof.* A line meets the open interior of a (closed convex) triangle iff the three affine
signs of the line at the triangle's three vertices contain both a strict $+$ and a
strict $-$. The **$(0,+,-)$ corner pattern** is a crossing: e.g. $y=0$ meets the
interior of the triangle bounded by $y=x$, $y=-x$, $x=1$, entering through the
multipoint $(0,0)$ (audited separately, `CrossingTheoremAudit` "chord endpoints at
concurrent triangle corners": $203{,}353$ such incidences checked exactly). The chord's
endpoints are $\ell\cap g$ for side-supporting lines $g$ of $T$ — finite vertices —
and they are distinct because $\ell$ does not support $T$ and $T$ is nondegenerate.
∎

**Lemma 6** (chord injectivity into unused gaps). *The open chord contains the relative
interior of at least one bounded elementary gap $e_t\subseteq\ell$; mapping each
crossing incidence to the leftmost such gap is injective, and no chord image is a side
image.*

*Proof.* Between the two distinct endpoint vertices of the chord lies at least one bounded
elementary gap, whose relative interior lies in the chord, hence in
$\operatorname{int}(T)$. Distinct triangles have disjoint open interiors; if two chords
output the same gap, the common relative interior lies in both open interiors, a
contradiction. If a chord output gap $e$ were a side image, then $e$'s relative
interior would lie inside a selected side $s\subseteq\partial T'$ *and* inside
$\operatorname{int}(T)$; choosing $x\in\operatorname{relint}(e)$, a small open disk
about $x$ lies inside $\operatorname{int}(T)$ (open set) while $T'$ contributes an
inward open half-disk there, so the open interiors meet — contradiction. This argument
is insensitive to boundary contact, overlapping boundaries, and shared segments. A side
on a line other than $\ell$ cannot contain a relative interval of an $\ell$-gap. Thus
chord images lie in $G_\ell\setminus U_\ell$, disjoint from side-gap images, and the
tokens are formal positions disjoint from gaps. ∎

**Theorem 4** (per-line inequality). *For each line $\ell$*
$$\sigma_\ell+\gamma_\ell\le m_\ell+2a_\ell.\tag{3}$$

*Proof.* Combine Lemma 4 (sides inject into gaps plus tokens) and Lemma 6
(crossings inject into gaps minus side-gaps). ∎

### 3.4 Counting and summation

**Lemma 7** (incidence identity). *$m_\ell+2a_\ell
= n-2-q_\ell-\sum_{p\,\text{on}\,\ell}(k_p-4)$, where the sum ranges over the
multipoints $p$ on $\ell$.*

*Proof.* The $n-1-q_\ell$ lines other than $\ell$ and not parallel to it each meet
$\ell$ once; grouping the $(k_p-1)$ such lines sharing the vertex $p$ gives
$n-1-q_\ell = r+\sum_{p\,\text{on}\,\ell}(k_p-2)$ (the $1+(k_p-2)$ split),
and so
$m_\ell+2a_\ell=r-1+2a_\ell
= n-2-q_\ell-\sum_{p\,\text{on}\,\ell}(k_p-2)+2a_\ell
= n-2-q_\ell-\sum_{p\,\text{on}\,\ell}(k_p-4)$. ∎

*Proof of Theorem 1.* Sum (3) over all $n$ lines, using Lemma 7. The left side
sums to $3|S|+C$. On the right, $\sum_\ell q_\ell=2Q$, and each multipoint $p$
lies on exactly $k_p$ lines, so the per-line inner sum contributes $k_p-4$ once
for each of those $k_p$ lines:
$\sum_\ell\sum_{p\,\text{on}\,\ell}(k_p-4)=\sum_p k_p(k_p-4)$. Hence
$3|S|+C\le n(n-2)-2Q-\sum_p k_p(k_p-4)$. ∎

### 3.5 The all-parallel exception is exactly the exception

**Lemma 8.** *Every line of an essential arrangement has a finite vertex. Conversely,
in the all-parallel arrangement $|S|=C=0$ while the right-hand side of (1) equals
$-n$.*

*Proof.* If $\ell$ had no finite vertex, every other line would be parallel to $\ell$,
making all lines parallel, contradicting essentiality. For the all-parallel class,
$Q=\binom{n}{2}$, there are no multipoints, and no nondegenerate triangle exists, so
the left side is $0$ while the right side is $n(n-2)-n(n-1)=-n<0$: (1) fails
for that single class, and it is the *only* class where the argument breaks (Lemma 1
uses $r\ge1$). Hence the essentiality hypothesis is minimal. ∎

### 3.6 A remark on terminology

Let $E_{\rm gap}=\sum_\ell m_\ell$ be the raw number of bounded elementary
arrangement segments. Grouping $k_p$ lines through $p$,
$E_{\rm gap}=n(n-2)-2Q-\sum_p k_p(k_p-2)$. The right-hand side of (1) is
$E_{\rm gap}+2\sum_p k_p$ — i.e. a *gap-plus-token* capacity, obtained from the raw
segment count by adding two formal slots for every line–multipoint incidence. This is
the numerically honest reading of the right-hand side of (1) (see
`CrossingTheoremAudit`, "independent_derivation_of_capacity"); the theorem itself is
unaffected.

---

## 4. Equality analysis and sharpness

### 4.1 When is the bound tight?

A family $S$ attains equality in (1) exactly when (a) every line $\ell$ attains
equality in (3); (b) no side image is ever displaced to a token except when the
gap-collision structure is exactly absorbed; and (c) the chord map and side map jointly
fill every bounded gap and every token on every line. Equality therefore requires, in
particular, that *every* bounded elementary segment of the arrangement is either covered by a
selected side or pierced by a selected triangle's interior, and that every multipoint
is used to its full two-token budget — a very rigid coincidental condition, which the
table below shows is met exactly twice among the known records.

### 4.2 The sharpness table

All entries were computed by **exact rational arithmetic** with the campaign's certified
routine `face_multiplicity_counts` (an exact separating-axis / sign-straddle oracle) on
*maximum interior-disjoint selections* of the published record arrangements, with
crossed triangles allowed. The arrangements at $n=9,11,13,15$ were machine-read
from Savchuk's LineOrder gallery (SVG combinatorics) and re-realized over
$\mathbb{Q}$ by this campaign's reconstruction pipeline; the $n=10$ row uses this
campaign's own exact-rational certificate (report.md §4.1). This table was
**independently reproduced** for this paper (`scratch/kobon/verify_sharpness_table.py`
→ `scratch/kobon/sharpness_table_verification.json`): every row passes
`engine.verify_selection` on the reconstructed exact arrangement and the
`face_multiplicity_counts` values below agree.

| $n$ | arrangement (source) | $|S|$ | $Q$ | multipoints | $C$ | $3|S|+C$ | RHS of (1) | slack |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 9 | Blanc-type $21$ (gallery `ibase_kobon_9.svg`) | 21 | 0 | none | 0 | **63** | 63 | **0** |
| 10 | this campaign's certificate (report.md §4.1) | 25 | 1 | none | 0 | 75 | 78 | 3 |
| 11 | record $32$ (gallery `kobon_11_32tri_lines.svg`) | 32 | 0 | none | 0 | 96 | 99 | 3 |
| 13 | record $47$ (gallery `kobon_13_m_sym_47tri_lines.svg`) | 47 | 0 | none | 0 | 141 | 143 | 2 |
| 15 | rotationally symmetric $65$ (gallery `kobon_15_5_rot_symmetry_lines.svg`) | 65 | 0 | none | 0 | **195** | 195 | **0** |

Arithmetic cross-checks: $9\cdot7=63=3\cdot21$; $15\cdot13=195=3\cdot65$;
$10\cdot8-2=78>75$; $11\cdot9=99>96$; $13\cdot11=143>141$.

**Corollary 4** (no constant strengthening). Since equality occurs at $n=9$ and
$n=15$, no constant can be subtracted from the right-hand side of (1) in general; the
bound is optimal as a universal integer inequality on the essential class.

**Gallery-source caveats.** The $n=9,11,13,15$ entries rest on geometry whose
*combinatorics* came from Savchuk's LineOrder gallery drawings. For such entries: (i)
the SVG is a drawing, not a certified proof; (ii) the reconstruction rationalizes
coordinates and re-verifies the *selection* exactly (disjoint interiors), so the
reported capacity numbers are exact *for the reconstructed arrangements*, which may differ
from the gallery intent in boundary cases (e.g. a near-parallel pair or a
near-concurrency could snap to exactness or to genericity); (iii) attribution metadata
("Blanc-type", "rotationally symmetric") is gallery-level and has not been independently
traced to a primary paper for every row; (iv) the values are *attained* quantities for
$|S|$, hence are unconditional lower bounds on $K_{\rm gen}(n)$, and the capacity
numbers are exact evaluations of (1) on those very configurations — the caveat concerns only
the reading "this is the *record* arrangement". The $n=10$ row carries no such caveat:
it is this campaign's exact certificate.

---

## 5. Comparison with prior bounds

### 5.1 What the literature contains

All prior upper bounds for Kobon-style counts live in one of two regimes, and neither
contains the quantities $C$ or $\sum_p k_p(k_p-4)$:

- **Tamura** — the bound $K(n)\le\lfloor n(n-2)/3\rfloor$ is attributed to Saburo
  Tamura via Wikipedia's and MathWorld's glosses of the Kobon problem. The statement is
  made for the broad problem, but the only locateable proofs are segment-charging
  arguments whose local steps are exactly the ones that fail in the broad convention
  (§5.3; `novelty_audit_capacity_theorem.md` row "Tamura"). No $Q$ term, no
  multipoint penalty, no $C$ term appears anywhere.
- **Clément–Bader** — *Tighter Upper Bound for the Number of Kobon Triangles*,
  unpublished ETH draft (21 Dec 2007; PDF cached at
  `oeis.org/A006066/a006066.pdf`). Theorem 1 gives the piecewise bound
  $K(n)\le\lfloor n(n-2)/3\rfloor-\mathbf 1_{n\bmod 6\in\{0,2\}}$ (with
  $(n-1)^2/3$ for $n\equiv1,4$). The proof is a charging argument whose two local
  indispensables are

  > (C–B Lemma 1 step 2) "If a line segment is the side of two triangles then the
  > corresponding line intersects at one of the two endpoints an existing point which belongs
  > to both triangles", and (step 3) "every intersection point $p$ with more than two
  > corresponding lines is part of at most two pairs of triangles that share a common side",

  both under a "perfect configuration" hypothesis. Neither statement survives in the broad
  convention; see §5.3 for exact witnesses.
- **BBL (Bartholdi–Blanc–Loisel)** — arXiv:0706.0723 (2007; see also Contemp.
  Math. 453, 2008). BBL bound the number $a_3(\mathcal{A})$ of triangular
  **faces** of **simple** arrangements, by definition ("two … intersect
  transversally into exactly one point, and this one does not belong to any other
  curve"), and prove in the affine case $a_3(\mathcal{A})\le \lfloor
  \frac{n(n-7/3)}{3}\rfloor$ for even $n$ (plus projective bounds). Because
  triangles are arrangement faces there, crossings, parallels and multipoints are
  **excluded by fiat**; consequently no simple-convention reference contains $C$ or the
  multipoint penalty — both are invisible in that regime. Their even-$n$ bound is
  *not* attained at $n=12,14$ pseudolines, as BBL themselves note.
- **Blanc (2011)**, *The best polynomial bounds for the number of triangles in a simple
  arrangement of $n$ pseudo-lines*, Geombinatorics **21** (2011) no. 1, 5–14
  (preprint arXiv:0801.2845, 2008) — same face convention. The introduction's
  hinge is that "a bounded segment
  may not delimit two different triangles, and the number of bounded segments is
  $n(n-2)$, [so] $a_3(\mathcal{A})\le n(n-2)/3$"; Corollary 2.4 gives
  $a_3(\mathcal{A})\le\lfloor n(n-5/2)/3\rfloor$ for even affine $n$ via an
  unused-segment parity argument (Prop. 2.3). Both hold for **simple** pseudoline
  arrangements only. The hinge sentence — one bounded segment serving at most one of the
  two triangles bordering it — is **false** for open-interior-disjoint selections in the
  broad convention: report.md §5's cautionary lemma found $256$ of $500$ random
  configurations with a segment double-covered from opposite sides and no multipoint
  endpoint (`scratch/kobon/audit_shard.py`; exact witness in §5.3(1) below).

The remaining located items are either constructions or different problems, and carry no
capacity inequality: Füredi–Palásti 1984 (PAMS 92, 561–566) builds arrangements
with $\sim n^2/3$ triangles; Roudneff 1996 (JCTB 66, 44–74) bounds the
*number of triangular faces* of **simple projective pseudolines**; Forge–Ramírez-Alfonsín
1998 (DCG 20, 155–161) gives a doubling construction only. The OEIS entry A006066
(Kobon sequence) tables constructions and bounds, and its companion *bound*
sequence is **A032765 = ⌊n(n+2)/3⌋**; no formula with $Q$ or $k_p$ penalties
appears there. Savchuk 2025 (arXiv:2507.07951) uses pseudoline-table encodings
that *can* represent parallels and multiple points and reports "no optimal solution
exists in the 11-line case" within that encoding, finding new optima at $n=23\,(161)$
and $n=27\,(225)$ — but its upper bounds are inherited from Tamura and Clément–Bader,
and it proves no degeneracy-refined inequality. Parpalak–Utkin 2026 (arXiv:2604.22035;
arXiv:2607.29236) study **bounded triangular faces in general position** only.

### 5.2 The scope gap

**report.md §2** asserts the located public tables record 25 and 38 as exact at $n=10,12$,
but without a proof chain that covers the broad convention: BBL/Blanc prove simple-face
statements (and BBL note their even-$n$ bound is not attained for 12 and 14
pseudolines); the Clément–Bader draft proves its charging in the perfect/nondegenerate
regime; Savchuk can represent parallels and multiple points in pseudoline tables but states no
exact theorem for $n\in\{10,12,14\}$; Parpalak–Utkin concerns bounded triangular faces
in general position. The novelty audit (`scratch/kobon/novelty_audit_capacity_theorem.md`)
checked each located source source-by-source and concluded:

1. **(Cap) is not subsumed.** No located published source proves any upper bound for
   open-interior-disjoint selections in the broad convention that is (a) refined by
   $Q$, $k_p$ and charged for crossed triangles by $C$. In the simple face convention
   the inequality degenerates to the classical $3T\le n(n-2)$ — a special case, not
   a subsumption. The Clément–Bader draft states only pointwise mod-6 integer bounds,
   via steps that are false in exactly the broad convention (witnesses below).
2. **(Cap) is not contradicted.** Nothing located asserts anything incompatible; the
   classical bounds agree with (1) in their common domain.
3. The novelty is *precisely* the $C$ term and the $\sum_p k_p(k_p-4)$ penalty —
   structurally invisible in the face convention ($C\equiv0$, all $k_p=2$).
4. The **$n(n-2)$ simple special case** is folklore; the claimable novelty is the
   broad-convention inequality **with** its two refinements and its proof.

The internally attempted weaker form $\sum_p(k_p-2)^2$ in place of
$\sum_p k_p(k_p-4)$ is **false** (report.md §6) — exact witnesses refute its local
charging; the surviving form is $(k_p)(k_p-4)$.

### 5.3 Why classical charging fails in the broad convention (exact witnesses)

The following three exact-rational configurations, each fully verified by the campaign's
exact routines (report.md §6), show that the two C–B local steps and the naive per-line
segment count cannot be transplanted:

**(1) A shared segment need not have a multiple endpoint** (kills C–B step 2 and its
"rescue"). $n=4$, lines $y=0$, $y=x/5$, $y=x-1$, $y=2x-4$; triangles
$\{0,1,2\}$ and $\{0,2,3\}$ are interior-disjoint and share the whole elementary
segment $(1,0)$–$(5/4,1/4)$ of $y=x-1$; both endpoints are simple crossings. No
parallels, no concurrency. This also kills the "at Tamura equality the excess $0$
forces full-side sharing" salvage: the line $y=x-1$ carries two selected sides yet is
fully tight.

**(2) A single multipoint can serve many shared-side pairs** (kills C–B step 3).
$n=6$, lines $y=0$, $y=x$, $y=-x$ (concurrent at the origin), $10y=10-x$, $x=-1$,
$20y=x-20$; the six triangles $\{0,1,3\},\{1,2,3\},\{2,0,4\},\{0,1,4\},
\{1,2,5\},\{2,0,5\}$ are pairwise interior-disjoint with **six** shared-side pairs
at the single triple point — in the face/perfect convention each multiple point serves at
most two.

**(3) The per-line segment count fails under concurrency.** $n=6$, lines $y=0$,
$y=x/5$, $x=1$, $10y=x-1$, $4y=x-2$, $3y=x-3$; the five pairwise-disjoint
triangles $\{0,1,2\},\{0,3,4\},\{0,4,5\},\{0,1,3\},\{0,2,4\}$ all have
a side on $y=0$, which carries only $4$ finite vertices / $3$ bounded elementary gaps:
$\sigma=5>n-2=4$. Hence the segment-count route cannot prove
$\lfloor n(n-2)/3\rfloor$ in the broad convention: $3T$ segments exceed $n(n-2)-2Q$
for a degenerate but perfectly legal arrangement. [Note: with multipoints the
correct per-line account is $\sigma_\ell+\gamma\le m_\ell+2a_\ell$, and (1) — this is
exactly what the theorem fixes.]

What *does* survive a straightforward proof (per-line lemma, no multipoint):
$\sigma_\ell\le k_\ell-1$ — "if no three lines are concurrent then
$3T\le n(n-2)-2Q$"; and the strict cotangent-descent argument. The generalized
identity (1) is what the present paper proves.

---

## 6. Consequences for the open targets

Specializing (2) gives uniform, audited *budget* constraints that the running SAT
campaign enforces as necessary clauses. Here $\Delta=n(n-2)-3T$.

### 6.1 $n=11$

At $T=33$: $\Delta=0$. Every hypothetical $33$-family must satisfy
$C+2Q+\sum_p k_p(k_p-4)\le0$. Since $k_p(k_p-4)=0$ at $k_p=4$ and $<0$ at
$k_p=3$, this forces $Q=0$ and no point of multiplicity $\ge5$; a triple point
contributes $-3$, so $C\le 3\,N_3$ where $N_3$ is the number of triple
points. In the **no-concurrency case the left side forces $C=0$: all 33 selected
triangles would have to be crossed-free — a *faces-only* 33-family** — which is the
(open) statement that no arrangement with parallels (and no triple points) beats the
face count at $n=11$. The known record is $32$ ($3\cdot32=96\le99$, slack $3$);
Savchuk (arXiv:2507.07951) reports a SAT case exhaustion within its pseudoline-table
encoding giving "no optimal solution exists in the 11-line case" there. The capacity
bound does not by itself decide $K(11)$.

### 6.2 $n=12$

At $T=39$: $\Delta=120-117=3$, so $C+2Q+\sum_p k_p(k_p-4)\le3$. Consequently:
$Q\le1$ without multipoints, $C\le3$ at $Q=0$ and $C\le1$ at $Q=1$; each
triple point adds $3$ units; $k_p=4$ points are budget-free; any $k_p\ge5$
consumes the whole budget. The record $K_{\rm gen}(12)\ge38$ (Kabanovitch-type
exact certificate with two triple points, `scratch/kobon/n12/n12_lower_certificate.json`)
attains $C=0$: with $Q=0$ and two $k_p=3$ vertices the RHS of (1) is
$120-2\cdot0-2\cdot(-3)=126$, so that family has slack $126-3\cdot38=12$; any
hypothetical $T=39$ family would need $C\le126-117=9$.

### 6.3 $n=14$

At $T=54$: $\Delta=168-162=6$, so $C+2Q+\sum_p k_p(k_p-4)\le6$. In the
*no-multipoint* branch ($\sum_p k_p(k_p-4)\ge0$, automatic when every
$k_p\ge4$) this forces $Q\le3$ strictly; with multipoints, $Q\ge4$ is not
excluded but requires compensation: each $k=3$ point *adds* $3$ units of
capacity ($k_p(k_p-4)=-3$), so e.g. $Q=4$ needs at least one triple point,
$Q=6$ at least two; $k=4$ points are budget-free and a $k=5$ point costs
$5$. The previous $53$ benchmark configurations live in the
no-crossing $Q\in\{2,3\}$ regime; the verified $54$ witnesses live at
$Q=0$ with shared-line triple points (Status below). The benchmark's two
**independently verified exact rational realizations** are the snapped
Bader-type realization
(`scratch/kobon/discoveries/n14_bader53_verified.json`, $Q=2$) and its
normalized sibling with parallel classes at positions $(6,7),(10,11),(12,13)$
(`scratch/kobon/discoveries/n14_best_exact53.json`, $Q=3$) — sitting at
$3\cdot53=159$ against the no-concurrency $Q=3$ ceiling $168-6=162=3\cdot54$: a
$54$-family inside that regime would saturate (1) with equality on every
line (the branch Theorem M rules out). The full $54$ case tree thus splits
into the no-concurrency Q-cubes ($q0,q1,q2,q3\text{paired},q3\text{triad}$)
and the multipoint-concurrency cubes — exactly the cube fleet under
construction. The campaign's exact certificate probes around the
snapped $Q=3$ realization are all UNSAT at radius $0..15$ (§8).

**Status (2026-08-20):** the broad-convention $53$-vs-$54$ question is settled
from below: Maiorana's exact-rational 14-line arrangement (source commit
`e47c7cfd`; independently re-verified here by three implementations including
`engine.verify_selection`; all 15 solutions re-verified by the replayable
sweep `scratch/kobon/n14/sweep_all_verify.py`, exit-0 transcript in
`sweep_all_verify.log`) witnesses $K_{\rm gen}(14)\ge54$ (§7.5). The witness
is non-simple ($N_3=2$ shared-line triples $(0,4,10)$, $(0,11,12)$, $Q=0$);
the remaining broad question is the **upper frontier $T=55$**
($\Delta=168-165=3$), carried by the all-degeneracy monolith. The
simple-face literature is separately scoped: neither the BBL upper bound nor
the non-simple Maiorana witness is used here to assert the face-convention
equality $K(14)=54$.

### 6.4 $n=18$ and $n=20$

At $T=94$ ($\Delta=288-282=6$) the $n=14$ tree applies verbatim: $Q\le3$,
$Q\ge4$ needs triple points, $k_p=4$ free, $k_p\ge5$ costly ($Q\le3$ with
$\sum k_p(k_p-4)\le6-2Q$). At $T=117$ ($\Delta=360-351=9$): $Q\le4$.

### 6.5 Reading the tree

None of these statements decides $K(n)$; (1) is a *necessary condition*. What it does
is replace an intractable search over all degeneracies by a finite, audited case tree
(the $n=14$ degeneracy lattice — $1{,}331$ admissible signatures $(Q,N_3,\dots,N_7)$
satisfying the capacity, bounded-face, and pair-account budgets — was enumerated exactly in
$0.4$ s, `scratch/kobon/n14/concurrency_cubes.md`). Each leaf is a cube of the
SAT cover.

---

## 7. The computational program: soundness, validation, and audits

### 7.1 Soundness architecture (report.md §3.4)

The engine reduces each target-$T$ question to SAT over a *combinatorial relaxation*
of line configurations, with variables for crossings, concurrencies, ordering along each
line, above/below relations, and selections, plus machine-verified lemmas (middle
crossing inequality, sidedness table, parallel linkage, disconnected-disjointness rules) —
every clause is **necessary**: it holds for every real configuration. Hence

- **UNSAT** at target $T$ proves $K_{\rm gen}(n)<T$ for *arbitrary* line
  configurations (no further hypothesis);
- **SAT** produces only a *candidate*: it is realized by fixed-slope LP straightening
  (order constraints become linear inequalities in intercepts; concurrencies become linear
  equalities solved exactly over $\mathbb{Q}$) and then **re-verified with exact rational
  arithmetic** — the verified configuration is the lower-bound certificate, independent
  of the solver.

The capacity inequality (1) enters as a pure necessary accelerator: the engine encodes per
line the exact-Tseitin family of
$\sigma_r+\gamma_r\le n-2-q_r-\sum_{p\in r}(k_p-4)$
(`add_face_bound(..., per_line=True)`), with a reduced no-concurrency variant
`add_no_concurrency_crossing_budget`; both are covered by the same necessity discipline
(17/17 TDD tests green).

### 7.2 Exact adversarial validation of the proof

The theorem's local injective steps were stress-tested, independently and exactly, at three
levels:

- **First audit** (`CrossingTheoremAudit`, verdict **VALID**, itemized by proof step):
  re-derives (1) independently; verifies gap-assignment chains, same-side injectivity
  (including nested bases and shared segments), the collision⇒multipoint step (with the
  explicit "no image collision but shared segment" counterexample bracketed), token
  injectivity (both endpoints multiple, adjacent collision gaps), chord injectivity under
  boundary contacts (including the $(0,+,-)$ corner), the incidence identity, and the
  edge cases $n=3$ (with a parallel pair; all-concurrent), a pencil plus parallels,
  all-concurrent, and shared selected edges — all with exact rational coordinates. It also
  notes the gap-plus-token reading of the RHS (§3.6).
- **Second audit** (`SuddenPtarmigan`, 58m37s): a dedicated adversarial harness
  (`scratch/kobon/audit_shard.py`) with independent oracles
  (`verify_selection`, `interiors_disjoint_h`, `face_multiplicity_counts`) tried to
  *break* the two weakest links on
  **$3{,}016$ arrangements** (integer, planted pencils with
  $k\in\{3,4,n-2,n-1\}$, planted concurrencies; $n\in\{5,\dots,8\}$),
  **$812{,}486$ valid families**, and
  **$5{,}900{,}863$ linewise checks**: **zero violations**; all $19$ violation flags at
  zero; $106{,}219$ tight lines saturating $\sigma+\gamma=m_\ell+2a_\ell$; the
  $(0,+,-)$ predicate agreed with the engine's incidence predicate on every one of the
  $812$k families; $58{,}174$ triangles enumerated; $527{,}667$ D∧A gap collisions
  (every one with a multipoint endpoint), $33{,}883$ with both endpoints multiple.
- **Novelty audit** (`novelty_audit_capacity_theorem.md`): no published source contains
  the $C$ term or the $\sum_p k_p(k_p-4)$ penalty; the two located near-misses
  (Tamura statement-only; C–B draft) are refuted by §5.3 witnesses in the broad
  convention. Verdict: **new — not subsumed, not contradicted**.

The proof itself is elementary and self-contained (this document); the audits are
independent *checks*, not prerequisites, of its validity.

### 7.3 Flagship certified application (context)

As the companion certified decision, $K_{\rm gen}(10)=25$: lower certificate (exact
rational §4.1 arrangement, $25$ triangles), upper certificate at target $26$ by an
exhaustive $11$-cube cover whose 9 complete DRAT proofs were byte-re-verified
(`drat-trim` `s VERIFIED`) and whose 2 truncated cubes were re-proven via 16
byte-verified subcubes (16 fresh DRAT proofs, all `s VERIFIED`; ≈ $39.8$ CPU-h
total). The boundary condition (1) is what made the $n=10$ structure tractable
(no-concurrency budget forcing $Q\le1$ at $T=26$).

### 7.4 In-flight targets (status, not claims)

- $n=11$ ($T=33$): faces-only consequence §6.1 in force; monolith under solving.
- $n=12$ ($T=39$): 15-cube cover, DRAT auto-verification pipeline running; the
  target remains **open** (discovery UNSAT is evidence, not certificate).
- $n=14$: the broad lower bound is **settled from below** — $K_{\rm gen}(14)\ge54$
  by the independently verified Maiorana witness (§6.3 status, §7.5); target
  $T=55$ (capacity $\Delta=3$) is the **open** upper frontier, with the
  all-degeneracy monolith and the Q-cubes at $T=54$ still completing as
  encoding-margin evidence.
- $n=18,20$: cubes built, solving.

Per the campaign's soundness rules, *no solver UNSAT is a theorem until DRAT-verified,
and no SAT model is a claim until exactly straightened and re-verified* — consequently
**$K(11)$, $K(12)$, $K(18)$, $K(20)$ are not asserted to be decided anywhere
in this document**, and the broad upper statement $K_{\rm gen}(14)\le54$ is
likewise not established (the $n=14$ lower bound is an exact-rational
certificate, not a solver SAT claim; its complement is the target $T=55$
frontier).

**Reserved verdict slots** (sections to be drafted when the corresponding runs
settle; nothing below is a claim).

**(a) $n=11$, $T=33$.** *Reserved.* A proof-grade run is in flight (watcher:
`N11ProofRelay`). This slot will record the certified outcome — verdict,
certificate chain (cube cover / DRAT verification / exact-rational model
realization, whichever applies), and its reading against the §6.1 budget
($C+2Q+\sum_p k_p(k_p-4)\le0$; faces-only in the no-concurrency case) — once
the run settles.

**(b) $n=14$ — final case-tree state.** *Reserved.* The broad lower bound
$K_{\rm gen}(14)\ge54$ is independently verified (§6.3, §7.5); the T=54
monolith and concurrency-cube fleet keep grinding as encoding-margin
evidence, and the **$T=55$ upper frontier** (all-degeneracy monolith) is
in flight. This slot will record the final case-tree state — per-cube
ledger, DRAT verification status, and the resulting position on
$K_{\rm gen}(14)$ vs $54$ under the §6.3 anatomy — once the fleets settle.

---

## 7.5 Results ledger and fill-in target

This section is the paper-facing status ledger. A row marked **open** is not an
upper-bound claim: under the campaign's soundness rule, a solver UNSAT becomes a
theorem only after an independently checked DRAT certificate, while a SAT model
becomes a lower-bound certificate only after exact-rational straightening and
selection verification. The budget recursion used throughout is

$$
M_n:=n(n-2),\qquad
\Delta_n(T):=M_n-3T,\qquad
\mathcal{B}_n(T):=C+2Q+\sum_{p}k_p(k_p-4),
\qquad \mathcal{B}_n(T)\le\Delta_n(T). \tag{4}
$$

Thus, in a no-concurrency, crossing-free branch, $2Q\le\Delta_n(T)$; this is
only a branch consequence, not a substitute for the full capacity theorem.
For rows with no current target, the dash is deliberate: no value is being
silently inferred from a public table.

### Results table

| $n$ | reference target / record | $\Delta_n$ | role and current status | pending verdict or certificate |
|---:|:---|---:|:---|:---|
| 9 | $T=21$ | $0$ | **SEALED control:** $K_{\rm gen}(9)=21$; the target-22 control proof is DRAT-verified and the 21-triangle model is exact-verified | none for the control; retain the archived recipe/proof hashes |
| 10 | $T=26$ | $2$ | **THEOREM:** $K_{\rm gen}(10)=25$; the exhaustive 11-cube cover is DRAT-certified, with the exact 25-triangle lower certificate | none |
| 11 | $T=33$ | $0$ | Lower certificate $K_{\rm gen}(11)\ge32$; the no-concurrency slice is discovery-UNSAT, while the 124-signature concurrency-inclusive monolith is still open | `[ ]` kissat UNSAT (or exact SAT realization) plus `drat-trim` verdict; only then fill $K_{\rm gen}(11)=32$ or the SAT alternative |
| 12 | $T=39$ | $3$ | Lower certificate $K_{\rm gen}(12)\ge38$; the 15-cube/monolith upper-bound front is open | `[ ]` complete the target-39 cube ledger and independently verify every DRAT proof |
| 13 | record $T=47$ | $2$ | Sharpness/control row: the exact reconstructed 47-triangle selection has $3T+C=141$ against capacity $143$; no global decision campaign is claimed | `[ ]` bibliography/provenance pin for the record; no $K_{\rm gen}(13)$ verdict claimed |
| 14 | $T=54$ | $6$ | **VERIFIED LOWER CERTIFICATE:** $K_{\rm gen}(14)\ge54$ from Maiorana's 14-line arrangement, independently verified by us on exact rational data (SHA-pinned sol1 `e47c7cf`, full sol1..sol15 sweep: 15/15 solutions; replayable sweep `scratch/kobon/n14/sweep_all_verify.py` + transcript `sweep_all_verify.log` (exit 0, persisted inputs in `scratch/kobon/n14/maiorana/` with manifest `SHA256SUMS.txt`); sol1 transcript `scratch/kobon/n14/maiorana_verify.log`; certificate `scratch/kobon/n14/Kgen14_ge54_certificate.md`) | This is a lower-bound certificate, not an upper-bound claim; target $T=55$ is the next broad-convention upper-bound object |
| 15 | record $T=65$ | $0$ | Sharpness/control row: the rotationally symmetric record attains equality in (1); this is not a new global exact-value claim | `[ ]` primary construction citation and any global-status statement |
| 16 | — | — | No target-specific front is registered in the current artifacts; the row is retained to prevent an unsupported interpolation between $n=15$ and $n=18$ | `[ ]` choose a target, lower certificate, and proof-relevant case tree |
| 17 | — | — | No target-specific front is registered in the current artifacts | `[ ]` choose a target, lower certificate, and proof-relevant case tree |
| 18 | $T=94$ | $6$ | Known exact lower certificate $K_{\rm gen}(18)\ge93$; the $Q=3$ no-concurrency branch follows the Theorem-M mechanism, while the built cube fleet is still solving | `[ ]` certify the remaining cubes (including q0/q1 and any residual concurrency branches) and close the target-94 ledger |
| 19 | — | — | No target-specific front is registered in the current artifacts | `[ ]` choose a target, lower certificate, and proof-relevant case tree |
| 20 | $T=117$ | $9$ | Known exact lower certificate $K_{\rm gen}(20)\ge116$; the budget tree has $Q\le4$ in the no-multipoint branch, and q0--q4 cube encodings are present; the target remains open | `[ ]` integrate the $Q=3$ obstruction (pair/triad branches), then settle the residual cubes and DRAT status |

### Budget table for the active targets

| front | arithmetic | immediate no-multipoint consequence | open-front reading |
|:---|:---|:---|:---|
| $n=11,T=33$ | $\Delta=99-99=0$ | $Q=0$ and $C=0$ when no triple point is present | the no-concurrency discovery cube is closed only at discovery tier; concurrency signatures remain |
| $n=12,T=39$ | $\Delta=120-117=3$ | $Q\le1$; at $Q=0$, $C\le3$, and at $Q=1$, $C\le1$ | capacity narrows the 15-cube cover but does not decide it |
| $n=14,T=54$ | $\Delta=168-162=6$ | $Q\le3$; $Q=3$ is the equality branch addressed by Theorem M | the exact-rational Maiorana witness establishes $K_{\rm gen}(14)\ge54$; target $T=55$ is the next broad-convention upper-bound front |
| $n=18,T=94$ | $\Delta=288-282=6$ | the same $Q\le3$ tree as $n=14$ | reuse the analytic $Q=3$ obstruction, but certify the larger residual encodings |
| $n=20,T=117$ | $\Delta=360-351=9$ | $Q\le4$ | the $Q=3$ pair/triad obstruction and the residual $Q=0,1,2,4$ fronts are separate checkboxes |

The entries above intentionally separate theorem recursion from verdict. In
particular, a negative or small $\Delta$ is not itself an upper bound in the
broad convention: triple points contribute the signed term
$k_p(k_p-4)$, and crossed triangles consume the $C$ budget.

## 7.6 Methods (draft)

### Encoding and scale

The SAT layer is a necessary-clause relaxation of straight-line arrangements.
The base model carries crossing, concurrency, order, parallel-orientation,
sidedness, betweenness, and selected-triangle variables; exact reified
crossing indicators, the crossing-refined per-line capacity clauses, the
bounded-face penalty, one audited $\sigma$ symmetry bit, and exact-$T$
selection are then added. The clauses are used only in the sound direction:
UNSAT is an upper-bound certificate after proof checking, whereas SAT is
decoded, fixed-slope straightened, and checked with exact rational
arithmetic.

For reproducibility, the four archived $n=9$ target-22 control
cube/recipe snapshots retained in the working bundle have these DIMACS
headers (variables / clauses):

| archived $n=9$ snapshot | variables | clauses |
|:---|---:|---:|
| `kobon_n9_t22_exact.cnf` | 5,384 | 254,101 |
| `kobon_n9_t22_exactface.cnf` | 68,091 | 384,729 |
| `kobon_n9_t22_canonface.cnf` | 8,672 | 262,045 |
| `kobon_n9_t22_canoncount.cnf` | 9,216 | 266,076 |

These are historical control encodings, not interchangeable proof objects.
The current positive target-21 rehearsal uses the explicit recipe
`build_model+TC+face_bound(per_line=True)+symbreak+exact` and has **62,517
variables / 379,298 clauses**; it yields the exact-verified 21-triangle lower
certificate. This distinction prevents the legacy 5,343-variable target-21
file from being mistaken for the current recipe.

The monolith and large-cube scale grows quickly. The byte-identified
$n=11,T=33$ monolith has **252,622 variables / 1,521,724 clauses**. The
current $n=12,T=39$ monolith log records **24,131 / 2,763,519**; the older
on-disk proof-CNF era is **60,194 / 1,925,748**, and must not be conflated
with the current file. The $n=14,T=54$ monolith parses as
**1,357,615 / 7,841,316** (the progress log rounds this to 1.36M / 7.8M).
For scale beyond the monoliths, the $n=18,T=94$ q0 cube is
**1,708,348 / 29,932,035**, and the $n=20,T=117$ q0 cube is
**3,455,871 / 59,195,137**; the other budget cubes differ mainly by their
fixed degeneracy units. These counts are DIMACS headers read from the
archived artifacts, not estimates inferred from file sizes.

### Solvers, seeds, and certificate discipline

Discovery and proof runs use **Kissat 4.0.4** and **CaDiCaL 1.9.5**. The
proof-grade recipe pins a recorded Kissat seed (the n=11 rerun is
`--unsat --seed=0`); discovery work also uses independent solver families and
seed racers where memory permits. Seed diversity is a robustness check, not a
replacement for a certificate, and the final per-cube seed/command manifest
remains a `[PENDING-VERIFICATION]` ledger item.

Every claimed UNSAT front must end in a DRAT proof whose `drat-trim` output is
`s VERIFIED`; a solver's exit code or a discovery log alone is not sufficient.
Every claimed SAT front must carry the decoded line arrangement, exact
straightening data, and exact `verify_selection` result. The $n=10$ cover and
the $n=9$ target-22 control satisfy this standard; the $n=11,12,14,18,20$
rows above are deliberately left as fill-in slots until the corresponding
proof chains settle.

## 7.7 Related Work (DRAFT; discrete geometry)

The Kobon problem belongs to discrete geometry at the intersection of line
arrangements, extremal triangle counts, and certified combinatorial search.
Savchuk (2025, [arXiv:2507.07951](https://arxiv.org/abs/2507.07951)) gives a
modern pseudoline-table/SAT encoding and closes the encoded $n\le11$ cases
(including the statement “no optimal solution exists in the 11-line case”);
the wording is encoding-relative and is not, by itself, a theorem for the
present broad straight-line convention. Blanc (2011, [arXiv:0801.2845](https://arxiv.org/abs/0801.2845);
[source PDF](https://algebra.dmi.unibas.ch/blanc/articles/bestbound.pdf)) studies
triangular faces of simple pseudoline arrangements and supplies the sharp
face-convention bounds, including the even affine $n=14$ value $53$. The
Clément--Bader (2007) ETH draft, [*Tighter Upper Bound for the Number of Kobon
Triangles*](https://www.sop.tik.ee.ethz.ch/publicationListFiles/cb2007a.pdf),
gives the familiar piecewise mod-$6$ upper-bound formula in its
perfect/nondegenerate regime; Bader's [14-line gallery entry](https://www.sop.tik.ee.ethz.ch/people/baderj/other.html) is
construction provenance, not a proof of the broad upper-bound frontier. The
campaign's Maiorana lower certificate is separate: the arrangement was
independently verified by us on exact rational data (SHA-pinned sol1
`e47c7cf`, full sol1..sol15 sweep: 15/15 solutions; replayable sweep
`scratch/kobon/n14/sweep_all_verify.py`, transcript
`scratch/kobon/n14/sweep_all_verify.log`; sol1 transcript
`scratch/kobon/n14/maiorana_verify.log`; certificate
`scratch/kobon/n14/Kgen14_ge54_certificate.md`), establishing
$K_{\rm gen}(14)\ge54$ under the broad convention. This does not assert the
standard simple-face value $K(14)=54$: the witness is non-simple, and the
face-convention literature is separately scoped.
[OEIS A006066](https://oeis.org/A006066)'s status convention is retained here: a row written with
“$\ge$” is a construction/lower-bound entry, while the adjacent number is a
reported bound or target, not an automatic exact value. In particular, the
currently approved row is `14 >= 53 54 [Bader]`; its approved construction
entry remains “$\ge53$ [Bader]”. The $54$ proposal in the update history is
unapproved OEIS history and is not used as an upper-bound claim here. The
source rows also give $18:\ \ge93,\ 94$ [Bader] and
$20:\ \ge116,\ 117$ [Wood]. Those sources frame the comparison made here:
the present capacity theorem allows parallels, multiple points, and crossed
selected triangles, and adds the explicit crossing charge $C$ and the signed
multipoint term $\sum_p k_p(k_p-4)$. Exact publication metadata, source-page
quotations, and the final bibliography ordering are **[PENDING: literature
pin]**; see `scratch/kobon/litrefs.md` and §5 for the current source audit.

### Open placeholders before submission

- [ ] Replace the n=11 discovery-only wording with the final DRAT-verified
  monolith verdict (or an exact SAT construction).
- [ ] Close and independently verify the n=12 target-39 fleet.
- [ ] Record the $n=14$, $T=55$ upper-bound cube/monolith ledger; the
  $T=54$ lower certificate is independently verified as documented above.
- [ ] Record the n=18 target-94 and n=20 target-117 proof chains, including
  the n=20 pair/triad obstruction citation and residual cube statuses.
- [ ] Decide whether the n=13, n=15 sharpness rows receive global $K_{\rm gen}$
  verdicts or remain construction/capacity controls; do not promote a record
  to an upper bound without a proof chain.
- [ ] Choose and document target fronts for $n=16,17,19$ if they remain in the
  paper's scope.
- [ ] Replace the four n=9 archive labels and the solver-seed sentence with
  the final frozen manifest if a later audit changes the control inventory.
- [ ] Insert the source-pinned Savchuk, Blanc, and Clément--Bader bibliography
  entries and verbatim status quotes; leave all unresolved attribution marked
  pending.

## 8. Limitations and future work

1. **Necessity only.** Condition (1) — and the SAT architecture around it — yields upper
   bounds and budget constraints; it does not construct or cut down to exact values.
   The exact Kobon numbers for $n\ge11$ remain open. In particular the $n=11$
   "33 forces faces-only" statement and the $n=14$ "$Q\le3$ in the
   no-multipoint branch" statement are *conditions any hypothetical extremal
   family in those respective branches must satisfy*; they do not establish the
   non-existence of such a family.
2. **No mod-6 / parity strengthening.** The classical face-convention refinements
   (the BBL/Blanc even-$n$ parity improvements and the Clément–Bader mod-6 piecewise
   formula) are of a qualitatively different nature than the degeneracy- and
   crossing-refined bound here, and are not consequences of (1). Obtaining a version
   of (1) that *also* captures parity of even-$n$ face counts is open.
3. **Straight lines only.** The proof is affine and uses slope/cotangent comparisons
   and straight-line intersection structure; it does **not** extend to general pseudoline
   arrangements, where the order can bend.
4. **Essentiality required.** The all-parallel class is the exact, unique exception
   (Lemma 8); no nondegenerate triangle exists there, and the inequality is false as
   literally stated (RHS $=-n$).
5. **Gallery limitations.** The sharpness rows for $n=9,11,13,15$ (§4.2) rest on
   exact reconstructions of drawings (gallery provenance); the values are exact *for the
   reconstructed arrangements* and are unconditional lower-bound capacities, but "record"
   attribution is gallery-level metadata.
6. **Prospective.** The $n=14$ lower-bound certificate $K_{\rm gen}(14)\ge54$
   (§6.3 Status, §7.5) shifts the nearest frontier: closing the target $T=55$
   all-degeneracy monolith (upper frontier, $\Delta=3$) and the certified
   closeout of $n=12$ are now the immediate next steps; the $Q=3$ placement
   cube is superseded workforce (Theorem M already closes that branch), and a
   strengthened per-line weighted token accounting, a projective-plane
   formulation, and $C$-weighted sharpenings remain the natural theoretical
   follow-ons.

---

## 9. Provenance and artifact index

- `math/kobon/report.md` — campaign report (theorem statement §7; certificate §4;
  soundness §3.4; failure-mode analysis §6; sharpness table §7).
- `math/kobon/engine.py` — audited engine; `face_multiplicity_counts`,
  `verify_selection`, capacity encodings; `AUDIT.md`, `README.md`.
- `history://CrossingBudgetMath` — construction of the proof (gap-and-token + chord).
- `history://CrossingTheoremAudit` — independent audit (VALID; itemized resolutions,
  airtight proof, edge cases, essentiality argument).
- `scratch/kobon/audit_shard.py` + `audit_shard_report.json` — second adversarial
  audit (exact $5{,}900{,}863$-check run, zero flags).
- `scratch/kobon/novelty_audit_capacity_theorem.md` — novelty audit with the
  citation-correction table (arXiv/identifier hygiene for BBL, Blanc, Roudneff,
  Clément–Bader, A032765).
- `scratch/kobon/verify_sharpness_table.py` +
  `sharpness_table_verification.json` — this paper's independent sharpness recomputation.
- `scratch/kobon/n14/concurrency_cubes.md` — $1{,}331$-signature degeneracy tree.
- `math/PROGRESS.md`, entry 2026-08-17 — settled/in-flight state of the campaign.

## References

[1] Fujimura, *The Tokyo Puzzles* (orig. 1978). [Kobon triangle problem; see
Wikipedia's list of unsolved problems in mathematics.]

[2] Report.md §2 sources as compiled by `scratch/kobon/novelty_audit_capacity_theorem.md`:
Tamura's bound as reported by Wikipedia/MathWorld (Kobon triangle problem;
KobonTriangle.html); the only locateable proofs are segment charging whose local steps fail in
the broad convention (§5.3).

[3] Clément–Bader, *Tighter Upper Bound for the Number of Kobon Triangles*, ETH draft,
21 Dec 2007 (unpublished); PDF https://www.sop.tik.ee.ethz.ch/publicationListFiles/cb2007a.pdf and cached at https://oeis.org/A006066/a006066.pdf.

[4] Bartholdi, Blanc, Loisel, arXiv:0706.0723 (2007); see also Contemp. Math. 453 (2008). Triangular faces of simple arrangements: affine $a_3(\mathcal{A})\le\lfloor n(n-7/3)/3\rfloor$ for even $n$; projective bounds.

[5] Blanc, *The best polynomial bounds for the number of triangles in a simple arrangement of $n$ pseudo-lines*, Geombinatorics **21** (2011) no. 1, 5–14; URL https://algebra.dmi.unibas.ch/blanc/articles/bestbound.pdf (preprint arXiv:0801.2845, 2008). Face-convention bound $a_3(\mathcal{A})\le n(n-2)/3$; even affine $a_3\le\lfloor n(n-5/2)/3\rfloor$; the segment hinge quoted in §5.

[6] Roudneff, *The maximum number of triangles in arrangements of pseudolines*, JCTB 66 (1996) 44–74.

[7] Forge, Ramírez-Alfonsín, *Straight line arrangements in the real projective plane*, DCG 20 (1998) 155–161.

[8] Füredi, Palásti, *Arrangements of lines with a large number of triangles*, PAMS 92 (1984) 561–566.

[9] OEIS A006066 (Kobon triangle sequence); companion bound sequence A032765 $=\lfloor n(n+2)/3\rfloor$.

[10] Savchuk, *Kobon triangles*, arXiv:2507.07951 (2025): constructions (new optima at $n=23\,(161)$, $n=27\,(225)$), "no optimal solution exists in the 11-line case" (within the encoding), upper bounds quoted from [2]/[3].

[11] Parpalak, Utkin, arXiv:2604.22035 and arXiv:2607.29236 (2026): bounded triangular faces in general position (constructions + non-degenerate enumeration; no degeneracy-refined inequalities).

[12] Alkauskas, *Triangle unions with maximal number of sides*, arXiv:2510.22584 (2025) — a different extremal problem (unions), noted as out of scope.

[13] Felsner, *Arrangements of pseudolines* (Handbook of Discrete and Computational Geometry, ch. 5), bibliography corroborating that the published "maximum triangles" theory is simple/face-convention only.

> **Citation care (from the novelty audit).** Two arXiv identifiers that have circulated
> in secondary sources are wrong and are *not* cited here: "Blanc 2011, arXiv:1012.1931"
> and "BBL arXiv:1705.07183" are particle-physics papers; the correct identifiers are
> Blanc arXiv:0801.2845 and BBL arXiv:0706.0723; the "Clément–Bader 0706.0726"
> identifier is likewise a physics paper and the true source is the ETH draft PDF [3].
> OEIS A008765 is unrelated; the companion *bound* sequence is A032765.

---

*End of paper. The numbered theorem, all numerical tables, all witness configurations and all
audit statistics in this document correspond 1:1 to the artifacts indexed in §9.*
