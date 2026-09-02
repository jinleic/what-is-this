# THEOREM H-DP1-THREEFACTOR — no genuinely three-dimensional circuit at $|S|=d+1$
## run `20260902T014115Z_25757c71_9c6065a4173c`

Agent `RsPe3dH2`, 2026-09-01/02 UTC. Original prereg commit/hash:
`1b904736706e7e248b723de05e945576f09d4fdf` /
`209c0f6a94e6fc5781b9d53853bcdb678b2e6ce94fe2e16891c7639194dbd442`.
Amended prereg commit/hash:
`0cc460f9100402c222964cd011afa80599df8b0e` /
`a34cf7a7622507815f0883eda549564e5b98e7fdd2f13457ec14515718c36b8c`.
Both byte-identical run copies are preserved and both hashes are recorded in
`provenance.txt`. The amendment narrowed one overbroad GRS-specific sentence;
§8 discloses it completely. The full unchanged T1–T5 family was rerun after
the amendment. No frozen directory was modified.

## 0. Setting and result

Let

$$H=H_1\otimes\cdots\otimes H_n,\qquad n\ge2,
\qquad d_i=\operatorname{spark}(H_i)\ge3,
\qquad d=\min_i d_i,$$

where every factor column is nonzero. The coordinate-$i$ projection of a grid
support $S$ is denoted $\pi_i(S)$; coordinate $i$ **varies** when
$|\pi_i(S)|>1$. Let $C_i(k)$ be the number of size-$k$ circuits of $H_i$.
For a pair $i<j$, let $N^{\rm nf}_{ij}(k)$ be the number of size-$k$ circuits
of $H_i\otimes H_j$ with both projections nontrivial. Thus `nf` means
nonfiber; this convention is required to avoid double counting.

**Theorem DP1-N.** Every size-$(d+1)$ circuit of $H$ varies in at most two
coordinates. If exactly one coordinate $i$ varies, the support is a transported
size-$(d+1)$ circuit of $H_i$. If exactly two coordinates $i,j$ vary, its
projection is a nonfiber size-$(d+1)$ circuit of $H_i\otimes H_j$ and every
other coordinate is fixed. Conversely every support of either form is a
circuit. Consequently

$$\boxed{
N_n(d+1)=
\sum_i\left(\prod_{k\ne i}s_k\right)C_i(d+1)
+\sum_{i<j}\left(\prod_{k\notin\{i,j\}}s_k\right)
N^{\rm nf}_{ij}(d+1).}
\tag{1}$$

There is therefore **no genuinely three-dimensional size-$(d+1)$ circuit**.
More sharply, the pair term can be nonzero only when $d=3$ and
$d_i=d_j=3$; then it is the sum of that pair's crossing profile-$(3,3)$ and
abstract all-distinct profile-$(4,4)$ populations. For two 2-row GRS factors,
the latter is the frozen bilinear/Möbius/cross-ratio family. No GRS structure
is assumed in the general theorem.

The candidate result had an honest missing premise: the previous two-factor
work had left repeated-mixed profiles open. Sections 1–3 close that gap before
using regrouping.

## 1. Minimal relations and the projection-floor lemma

Let $A\otimes B$ be a two-factor product, with nonzero columns and finite
sparks $d_A,d_B\ge3$. Put $d=\min(d_A,d_B)$ and let $S$ be a size-$(d+1)$
circuit. Write

$$U=\pi_A(S),\quad V=\pi_B(S),\quad p=|U|,\quad q=|V|.$$

Because $S$ is minimally dependent, its rank is exactly $|S|-1=d$ and its
relation space is one-dimensional. In particular a relation

$$\sum_{(u,v)\in S}\gamma_{uv}
  a_u\otimes b_v=0
\tag{2}$$

has every $\gamma_{uv}\ne0$: a zero coefficient would give a dependence on a
proper subset.

**Lemma 1 (projection floor).** Either $p=1$ and $S$ is a $B$-fiber, or
$p\ge d_A$. Symmetrically, either $q=1$ and $S$ is an $A$-fiber, or
$q\ge d_B$. Hence every nonfiber circuit has

$$p\ge d_A,\qquad q\ge d_B.
\tag{3}$$

*Proof.* Group (2) by its $A$ index:

$$\sum_{u\in U}a_u\otimes w_u=0,
\qquad
w_u=\sum_{v:(u,v)\in S}\gamma_{uv}b_v.
\tag{4}$$

Suppose $p<d_A$. By the definition of spark, the distinct columns
$\{a_u:u\in U\}$ are independent. Dual functionals isolating these columns in
their span, applied to (4), force $w_u=0$ for every occupied row. Since every
coefficient in (2) is nonzero, each row is then a nonzero relation among its
distinct $B$ columns and must contain at least $d_B$ cells. If $p\ge2$, this
gives

$$|S|=\sum_{u\in U}\deg(u)\ge2d_B\ge2d>d+1,$$

which is impossible. The only remaining case is $p=1$, a fiber. Tensoring by
the fixed nonzero $A$ column is injective (choose a dual functional taking
that column to one), so the projected $B$ set is itself a size-$(d+1)$ circuit.
The symmetric argument proves the other half. $\square$

This lemma already eliminates the putative profiles $(2,d)$ and their
transposes; it does not rely on a kernel-summand decomposition or a choice of
such a decomposition.

## 2. Rank-sum inequality, including repeated profiles

Enumerate the $m=d+1$ cells of $S$ as edges $e=1,\dots,m$. Let
$A_E=[a_{u_e}]$ and $B_E=[b_{v_e}]$, retaining repeated columns when an index
occurs on more than one edge. Equation (2) is exactly

$$A_E\operatorname{diag}(\gamma)B_E^{\mathsf T}=0.
\tag{5}$$

The diagonal is invertible. Every row of
$A_E\operatorname{diag}(\gamma)$ is therefore a relation on the edge columns
of $B_E$. Its row space, of dimension $\rho_A:=\operatorname{rank}A_E$, lies
inside a relation space of dimension $m-\rho_B$, where
$\rho_B:=\operatorname{rank}B_E$. Thus

$$\rho_A+\rho_B\le m=d+1.
\tag{6}$$

This argument does **not** require all indices to be distinct; it is precisely
the missing repeated-profile tool.

For a nonfiber circuit, Lemma 1 gives at least $d_A$ distinct $A$ columns and
at least $d_B$ distinct $B$ columns. Every $d_A-1$ distinct $A$ columns are
independent, and similarly on $B$, so

$$\rho_A\ge d_A-1,
\qquad \rho_B\ge d_B-1.
\tag{7}$$

Combining (6)–(7),

$$d_A+d_B\le d+3.
\tag{8}$$

But $d_A,d_B\ge d\ge3$. Hence $2d\le d+3$, so $d\le3$. The standing floor
forces $d=3$, and (8) then forces

$$d_A=d_B=3.
\tag{9}$$

**Consequence.** If $d\ge4$, every size-$(d+1)$ two-factor circuit is a
fiber. If $d=3$ but one factor has spark greater than three, every size-four
circuit is again a fiber. All nonfiber residue is confined to two spark-three
factors at four cells.

## 3. Exhausting the four residual profiles

Under (9), Lemma 1 and $|S|=4$ leave exactly

$$(p,q)\in\{(3,3),(3,4),(4,3),(4,4)\}.$$

- **$(3,3)$.** Here $|S|=d_A+d_B-2$. The frozen Theorem X converse,
  including its H2 closure, says exactly that $S$ is a crossing support built
  from one size-three circuit of each factor.
- **$(4,4)$.** This is exactly the abstract all-distinct two-factor channel.
  Nothing further is needed for the general theorem or count (1). When both
  factors are 2-row GRS matrices, the frozen successor theorems identify it
  with a nondegenerate bilinear/Möbius graph and evaluate it by the PGL sum.
- **$(3,4)$ is impossible.** Because there are four edges and four distinct
  $B$ vertices, every $B$ vertex has degree one. The three $A$ degrees are
  therefore $(2,1,1)$. Inequality (6), together with spark three on both
  sides, forces $\rho_A=\rho_B=2$. The three selected $A$ columns consequently
  have a one-dimensional kernel spanned by a vector
  $r=(r_0,r_1,r_2)$ with all entries nonzero (any two columns are independent).

  Group (2) into the three vectors $w_0,w_1,w_2$ as in (4). Exactness of
  tensoring over a field gives

  $$\ker(A_U\otimes I)=\ker(A_U)\otimes F^{r_B};$$

  hence there is one vector $z$ such that $w_i=r_i z$ for all three rows. Let
  rows 1 and 2 be the singleton rows, containing distinct columns
  $b_{v_1},b_{v_2}$. Then

  $$\gamma_1 b_{v_1}=r_1z,
  \qquad \gamma_2 b_{v_2}=r_2z.$$

  All four displayed scalars are nonzero, so $z\ne0$ and the two distinct
  $B$ columns are proportional. This contradicts spark$(B)=3$, which makes
  every pair independent.
- **$(4,3)$ is impossible** by the symmetric argument.

We have now proved the complete two-factor classification that the candidate
regrouping needed:

> A size-$(d+1)$ circuit for two factors of spark at least three is either a
> transported factor circuit (fiber), or $d_A=d_B=d=3$ and it is a crossing
> $(3,3)$ support or an all-distinct $(4,4)$ support. There are no other
> profiles.

This closes the predecessor's registered “repeated-mixed” open item in the
entire size-$(d+1)$ window, not merely at the swept configurations.

## 4. Induction on the number of factors

We prove DP1-N by induction on $n$. The two-factor case is §§1–3. For the
inductive step regroup

$$A=H_1,
\qquad B=\bigotimes_{j=2}^n H_j,
\qquad d_B=\min_{j\ge2}d_j.$$

All columns of $B$ are nonzero. Frozen T-DGE gives

$$\operatorname{spark}(B)=d_B$$

and says every size-$d_B$ circuit of $B$ is an axis fiber in exactly one
constituent factor. Apply the complete two-factor classification to
$A\otimes B$.

### 4.1 Fiber cases

If the composite $B$ index is fixed, the support transports a size-$(d+1)$
circuit of $H_1$ and varies in only coordinate 1.

If the $H_1$ index is fixed, the support projects to a size-$(d+1)$ circuit of
$B$. There are three exhaustive cases, which are necessary when the regrouped
sparks are unequal:

1. $d_B=d$: its size is $d_B+1$, so the induction hypothesis applied to the
   $n-1$ factors inside $B$ gives at most two varying coordinates;
2. $d_B=d+1$: its size equals spark$(B)$, so T-DGE makes it an axis fiber and
   it varies in one coordinate;
3. $d_B>d+1$: its size is below spark$(B)$, so this case is impossible.

Thus a fiber of the regrouped product never creates three-coordinate support.

### 4.2 Crossing case

A regrouped crossing can occur only when
$d=d_1=d_B=3$. Its $B$ projection is a size-three circuit. T-DGE makes that
triple an axis fiber along one coordinate $j\ge2$, with every other coordinate
of $B$ fixed. Therefore the global crossing varies only in coordinates
$1,j$.

### 4.3 All-distinct case

A regrouped all-distinct support likewise can occur only at
$d=d_1=d_B=3$. Let its four distinct $B$ columns have rank $\rho_B$. The
rank-sum inequality (6) applies to this very relation. The four distinct
$H_1$ columns have rank at least two because spark$(H_1)=3$, so

$$\rho_B\le4-2=2.$$

Every three-subset of the four $B$ columns is therefore dependent. Since
spark$(B)=3$, every pair is independent; each triple is consequently a
genuine size-three circuit. By T-DGE each triple is an axis fiber in $B$.
Take two triples sharing two columns, for example the triples on
$\{b_1,b_2,b_3\}$ and $\{b_1,b_2,b_4\}$. Two distinct product-grid points
determine their axis line: they differ in one coordinate, and any axis fiber
containing both must vary in that same coordinate and fix every other
coordinate to their common values. The two triples therefore lie on the same
axis line, and all four $B$ columns lie on it. The global all-distinct support
varies only in coordinate 1 and that one constituent coordinate $j$ of $B$.

In the crossing and all-distinct cases, tensoring the induced
$H_1\otimes H_j$ columns by the fixed nonzero columns of all other factors is
an injective linear transport: apply dual functionals taking each fixed column
to one. Dependence and every proper-subset independence are preserved. Thus
the two-coordinate projection is a genuine nonfiber circuit of
$H_1\otimes H_j$, not merely a support with the right profile.

The fiber, crossing, and all-distinct cases exhaust the regrouped two-factor
classification. The induction is complete. $\square$

## 5. Exact count and the required `nonfiber` convention

By §4 every circuit has a unique varying-coordinate set of cardinality one or
two.

- For varying set $\{i\}$, choose one of the $C_i(d+1)$ factor circuits and
  independently fix one index in every other factor. This gives
  $(\prod_{k\ne i}s_k)C_i(d+1)$ supports.
- For varying set $\{i,j\}$, choose one of the
  $N^{\rm nf}_{ij}(d+1)$ two-factor nonfiber circuits and independently fix
  every other coordinate. This gives
  $(\prod_{k\notin\{i,j\}}s_k)N^{\rm nf}_{ij}(d+1)$ supports.

Different varying-coordinate sets cannot describe the same support. Within a
fixed class, the varying projection and fixed indices recover the choices
uniquely. Conversely, tensoring either chosen circuit by fixed nonzero columns
is injective and preserves minimal dependence, so every constructed support is
a circuit. The classes are disjoint and exhaustive, proving (1).

If one instead defines $N_2$ to include pair fibers and then adds the
single-factor sum, (1) is false: a one-axis fiber appears through every pair
containing that axis. The exact formula requires $N^{\rm nf}_{ij}$. The T3
plant measured this failure explicitly: the correct 624-support family became
656 raw emissions, with 16 single-axis supports each appearing three times and
32 excess emissions.

A useful simplification follows from §§2–3:

$$N^{\rm nf}_{ij}(d+1)=0$$

unless global $d=3$ and $d_i=d_j=3$. For a surviving pair,

$$N^{\rm nf}_{ij}(4)=N^{\rm crossing}_{ij}+N^{\rm all}_{ij},$$

where the first term is the frozen crossing count and the second is the actual
all-distinct pair count (the PGL sum in the 2-row GRS specialization).

## 6. Exhaustive exact evidence

The standalone instrument used exact Python integers modulo 7, actual nested
row-major Kronecker columns, and exact Gauss–Jordan ranks. A support was called
a circuit only when its full column set had rank $m-1$ and every deletion had
rank $m-1$. It independently constructed the theorem family from actual factor
circuits plus actual **nonfiber** two-factor circuit censuses, fixed all other
indices, and compared full support sets. All four measured/constructed sets
were exactly equal; constructed multiplicity was one everywhere.

### T1 — $V_3^{\otimes3}$

Each $V_3$ is the $2\times3$ Vandermonde factor on points 1,2,3 over GF(7).

- Exhausting all $\binom{27}{3}=2925$ triples found exactly 27 size-three
  circuits, all fibers: profiles $(3,1,1),(1,3,1),(1,1,3)$ had 9 each.
- Exhausting all $\binom{27}{4}=17550$ four-subsets found exactly 81 circuits:
  profiles $(3,3,1),(3,1,3),(1,3,3)$ had 27 each.
- Each of the three two-factor references had exactly nine nonfiber crossings;
  fixing the third index gave $3\times3\times9=81$.
- No circuit varied in all three coordinates, and the 81-support sets were
  exactly equal.

### T2 — $V_4^{\otimes3}$

Each $V_4$ is the $2\times4$ Vandermonde factor on points 1,2,3,4 over GF(7).
The instrument exhausted all $\binom{64}{4}=635376$ supports.

| profile | count |
|---|---:|
| $(1,3,3),(3,1,3),(3,3,1)$ | 576 each |
| $(1,4,4),(4,1,4),(4,4,1)$ | 32 each |

Total: 1824. Each two-factor reference was independently 152 = 144 crossing
+ 8 all-distinct; fixing the third coordinate gives
$3\times4\times152=1824$. Zero circuit varied in all three coordinates, and
the full 1824-support census equaled the constructed set exactly.

### T3 — $M_{34}\otimes V_4\otimes V_4$

Here $M_{34}$ has columns $(1,x,x^2)^\mathsf T$, $x=1,2,3,4$, measured spark
4, and its unique four-subset is a genuine size-four circuit. This explicitly
tests a factor carrying a $(d+1)$-circuit and unequal regrouped sparks. All
635376 supports were exhausted.

| source/profile | count |
|---|---:|
| $M_{34}$ circuit, $(4,1,1)$ | $4\cdot4=16$ |
| $V_4\otimes V_4$ crossings, $(1,3,3)$ | $4\cdot144=576$ |
| $V_4\otimes V_4$ all-distinct, $(1,4,4)$ | $4\cdot8=32$ |
| **total** | **624** |

The two unequal-spark pair references each had four fiber circuits and **zero
nonfiber** circuits; the tied pair had 152 nonfiber circuits. The measured and
constructed 624-support sets were equal, with zero genuinely 3-D support.

### T4 — $M_{34}\otimes V_3\otimes V_4$

All $\binom{48}{4}=194580$ supports were exhausted. Exactly 156 circuits
occurred: 12 factor fibers of profile $(4,1,1)$ and
$4\times36=144$ $V_3\otimes V_4$ crossings of profile $(1,3,3)$. Both
unequal-spark pairs had zero nonfiber circuit. There was no other profile, no
genuinely 3-D support, and exact measured/constructed set equality held.

Together T3 and T4 test both branches that the equal-spark MDS anchors do not:
a real single-factor $d+1$ summand, an unequal composite spark, unequal factor
sizes, and the separation of pair fibers from $N^{\rm nf}$.

## 7. Controls in both directions

- **Known-true ACCEPT.** A concrete $V_3\otimes V_3$ crossing, embedded at one
  fixed third coordinate, had actual three-factor rank 3, all deletion triples
  independent, and was emitted by the constructor.
- **Non-tensor REJECT.** With the same labels, the first column was replaced by
  itself plus the true product column having the same first two indices and a
  distinct third index. Actual rank rose to 4 while the label constructor still
  accepted; the rank path rejected it.
- **Genuine-3D REJECT.** The four $V_4^{\otimes3}$ diagonal cells
  $(t,t,t)$ have profile $(4,4,4)$ but actual rank 4. They were absent from both
  measured circuits and constructed supports.
- **Duplicate REJECT.** Cloning a $V_4$ factor column dropped measured spark
  from 3 to 2 and exposed exactly the planted dependent pair $(0,4)$.
- **Three-factor Kronecker guard.** Nested product coordinates had dimension 8
  and matched the independent triple-loop coordinates exactly. Block
  concatenation had dimension 6; dropping the third factor had dimension 4;
  both were rejected.
- **Count-semantics REJECT.** Including pair fibers in $N_2$ produced 656 raw
  emissions rather than 624, duplicate excess 32, with 16 supports of
  multiplicity 3. This catches precisely the ambiguous version of the
  candidate formula.

The post-amendment run passed 78 aggregate assertions with zero failures. CPU
was 24.282/600 seconds and wall time 26.389/900 seconds at nice 15, all four
thread caps 1, exact arithmetic only, and bytecode disabled.

## 8. Disclosed prereg defect and resolution

After the first complete control pass, a proof audit found that the original
P1 sentence invoked the frozen cross-ratio/PGL result for every arbitrary
profile-$(4,4)$ factor pair. That result assumes two 2-row GRS factors, whereas
DP1-N assumes only nonzero columns and sparks at least three. This was a scope
error in one proof-obligation sentence, not a failed numeric claim.

The correction is structural and exact: for arbitrary factors, $(4,4)$ is
retained as the abstract all-distinct two-factor population inside
$N^{\rm nf}_{ij}$; the Möbius/PGL normal form is asserted only for the GRS
anchors. Sections 3–5 show that no normal form for that population is needed:
regrouping forces its composite-side columns onto one T-DGE axis, and count (1)
uses the actual pair population.

The original prereg and first-pass `controls_results.json` are preserved as
`pre_statement_original_preserved.md` and
`controls_results_pre_prereg_amendment_preserved.json`. The amended prereg was
committed, copied byte-identically to `pre_statement.md`, and provenance names
both commits and both hashes. `defect_log.json` records cause, impact, and
resolution. No target, pin, instrument, theorem statement, or promotion rule
was weakened. The complete unchanged T1–T5 suite was rerun after amendment and
passed as reported in §§6–7.

## 9. Scope and verdict

**Certified:** for every finite tensor product of factors with nonzero columns
and $d_i\ge3$, the complete size-$(d+1)$ support theorem DP1-N, including the
absence of genuinely three-coordinate circuits; the two-factor repeated-mixed
exclusion; and exact disjoint count (1). The analytic proof is field-independent
within the standing frozen T-DGE/Theorem X hypotheses; the numerical controls
are exact over GF(7).

**Not certified here:** spark-2 factors; any size larger than $d+1$; a closed
formula for arbitrary all-distinct pair populations beyond the 2-row GRS PGL
case; or a general circuit classification at other cardinalities for three or
more factors.

The candidate's honest completeness gap is closed in §§1–3, the induction's
unequal-spark cases are explicit, all exhaustive support sets equal the
constructed families, every supplied and added anchor matches, and every plant
fires. The resolved prereg scope amendment changes no theorem or evidence and
is fully provenance-bound. The unique supported verdict is
**FROZEN-CERTIFIED**.
