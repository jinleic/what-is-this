# THEOREM H-DP2-SIZE5 — complete size-five classification for two 2-row GRS factors
## run `20260902T021721Z_7629001e_215016bf40c3`

Agent `RsPe3dH2`, 2026-09-01/02 UTC. Prereg:
`prereg/H_DP2_SIZE5_PREREG_2026-09-01.md`, path-scoped commit
`515a3463110632f2d5d794983d4cbd4ff79dd930`, sha256
`dea07e216113a1cc90cb1d5a69e824360efc803267f802d6e1588eb3babdb491`,
byte-identical run copy `pre_statement.md`. No parameter-dependent computation
preceded that commit, no prereg amendment occurred, and no frozen directory or
root ledger was modified.

## 0. Setting and theorem

Let $A,B$ be 2-row GRS/Vandermonde parity-check factors over a field $F$, with
columns

$$a_x=(1,x)^{\mathsf T},\quad x\in X,
\qquad b_y=(1,y)^{\mathsf T},\quad y\in Y,$$

where the evaluation points in each set are distinct. Thus
$d_A=d_B=3$, and the product columns are

$$h(x,y)=a_x\otimes b_y=(1,y,x,xy)^{\mathsf T}\in F^4.$$

Write the profile of a grid support $S$ as
$(|\pi_XS|,|\pi_YS|)$. A path $P_k$ has $k$ vertices and $k-1$ edges.

**Theorem SIZE5.** A five-cell support $S\subseteq X\times Y$ is a circuit iff
all of the following hold:

1. every row and column degree in its bipartite incidence graph is at most two;
2. it contains no four-edge crossing support;
3. every all-distinct four-edge subset $Q\subset S$ has
   $\det[h(q)]_{q\in Q}\ne0$.

The possible circuit templates and their tests are exactly:

| profile | circuit graph | remaining test |
|---|---|---|
| $(3,3)$ | $C_4+K_2$ | automatic |
| $(3,4)$ | $P_5+K_2$ | automatic |
| $(3,5)$ | two row-centered $P_3$'s $+K_2$ | automatic |
| $(4,3),(5,3)$ | transposes of the preceding two | automatic |
| $(4,4)$ | $P_4+2K_2$ | one all-distinct 4-minor nonzero |
| $(4,5)$ | one row-centered $P_3+3K_2$ | two all-distinct 4-minors nonzero |
| $(5,4)$ | transpose | two all-distinct 4-minors nonzero |
| $(5,5)$ | $5K_2$ | all five 4-minors nonzero |

The six oriented tight cancellation channels $(1,2)$ and $(2,1)$ give five
distinct field-free profiles (the $(3,3)$ support is common to both
orientations). Their counts for $|X|=|Y|=n$ are

$$\boxed{N_{33}=9\binom n3^2,}$$

$$\boxed{N_{34}=N_{43}=18n\binom n3\binom{n-1}3
=72\binom n3\binom n4,}$$

$$\boxed{N_{35}=N_{53}=18n\binom n3\binom{n-1}4.}$$

For arbitrary $X,Y$, let

$$M_4=24\binom{|X|}4\binom{|Y|}4$$

be the number of all-distinct four-matchings and let
$N_{\rm all}^{(4)}(X,Y)$ be the frozen size-four all-distinct circuit count.
Then the remaining profile with a closed reduction is

$$\boxed{N_{44}^{(5)}=
12\bigl(M_4-N_{\rm all}^{(4)}(X,Y)\bigr).}$$

Profiles $(4,5),(5,4),(5,5)$ have the exact finite determinant sums in §7;
they are genuinely field-dependent and no field-free polynomial is claimed.

## 1. The exact support-cardinality law

Use the certified kernel identity to write any relation matrix as

$$\Gamma=C_1+C_2,$$

where every nonzero column of $C_1$ belongs to $\ker A$ and every nonzero row
of $C_2$ belongs to $\ker B$. Put

$$P=\operatorname{supp}C_1,
\quad Q=\operatorname{supp}C_2,
\quad T_1=|P|,
\quad T_2=|Q|,
\quad r=|P\cap Q|,$$

and let $c$ be the number of overlap cells where the two nonzero entries cancel
exactly. Outside $P\cap Q$, a nonzero entry of one summand survives. Inside the
overlap, a cell survives once if it does not cancel and zero times if it does.
Therefore

$$\boxed{|\operatorname{supp}\Gamma|=T_1+T_2-r-c.}\tag{1}$$

This corrects the candidate law's scope. If there are exactly $k_1$ and $k_2$
minimum-support factor relations, so
$T_1=k_1d_A$, $T_2=k_2d_B$, **and every overlap cancels**, so $r=c$, then and
only then (1) specializes to

$$|S|=k_1d_A+k_2d_B-2c.\tag{2}$$

Formula (2) is not the unrestricted law. The profile-$(4,4)$ circuits below
provide a load-bearing counterexample to that misuse:
$(T_1,T_2,r,c)=(6,6,4,3)$ gives $6+6-4-3=5$ by (1), whereas blindly using
$6+6-2c$ gives 6.

### Tight irreducible size-five shapes

Assume now $d_A=d_B=3$, every active relation has minimum support three, every
overlap cancels, and the decomposition is irreducible: each active row and
column contributes at least one survivor. Each active minimum relation can
then contain at most two overlap cells, so

$$c\le2k_1,\qquad c\le2k_2.\tag{3}$$

Equation (2) at size five reads

$$c=\frac{3(k_1+k_2)-5}{2}.\tag{4}$$

Suppose without loss $k_1\le k_2$. From (3)–(4),

$$3k_2-5\le k_1.$$

Since $k_2\ge k_1$, this implies $2k_1\le5$, so $k_1\le2$. If $k_1=1$,
parity in (4) and the inequality force $k_2=2$, giving $c=2$. If $k_1=2$,
the inequality forces $k_2=2$, but then (4) is nonintegral. Thus the only
irreducible tight/full-cancel solutions are

$$(k_1,k_2,c)=(1,2,2),\qquad(2,1,2).\tag{5}$$

Formal larger integer solutions necessarily have an active minimum relation
with no survivor; they are padded descriptions, not irreducible circuit
shapes. The other profiles must, and do, use a surviving overlap or a factor
relation of support greater than three.

## 2. Why the graph/minor criterion is complete

Five columns in the four-dimensional product space are always dependent. They
form a circuit exactly when every four-subset is independent: that condition
forces rank at least four, hence rank exactly four, and independence is
hereditary to all smaller subsets.

A row or column of degree at least three contains a size-three axis fiber.
Because every factor triple is a circuit, such a proper subset is dependent;
so a five-circuit has maximum degree two on both sides. Consequently both
profile coordinates are at least $\lceil5/2\rceil=3$.

Once degree-three fibers are excluded, any dependent four-subset is itself a
size-four circuit. The frozen complete size-four theorem says the only such
circuits are:

- profile $(3,3)$ crossings;
- profile $(4,4)$ all-distinct supports with zero bilinear determinant.

The formerly open repeated profiles $(3,4),(4,3)$ have been analytically
excluded. Therefore a five-set has no dependent proper subset exactly when it
has maximum degree two, no crossing four-subset, and no zero all-distinct
four-minor. This proves the stated criterion.

## 3. Exhausting the bipartite graph templates

A graph of maximum degree two is a disjoint union of paths and even cycles.
There are five edges; all listed projection vertices are occupied.

### Profiles with a 3-coordinate

- **$(3,3)$.** There are six vertices and five edges. A connected component is
  the tree $P_6$. Otherwise the only possibility is $C_4+K_2$. Deleting the
  middle edge of $P_6$ leaves two oppositely centered $P_3$ components, exactly
  a crossing; hence $P_6$ is rejected. $C_4+K_2$ contains no crossing and no
  all-distinct four-subset, so it is always a circuit.
- **$(3,4)$.** There are seven vertices and five edges. A $C_4$ would leave one
  edge unable to occupy the remaining one row and two columns, so the graph is
  a two-component forest. The edge partitions are $(4,1)$ and $(3,2)$:
  $P_5+K_2$ or $P_4+P_3$. In the latter, deleting the appropriate end edge of
  $P_4$ produces two oppositely centered $P_3$ components, a crossing. The
  former has neither forbidden minor and is always a circuit.
- **$(3,5)$.** Every column has degree one, while the row degrees are
  $(2,2,1)$. The unique graph is two row-centered $P_3$ components plus
  $K_2$. Its two $P_3$ centers lie on the same side, so they are not a crossing;
  there is no all-distinct four-subset because only three rows occur. It is
  always a circuit.

Transposition gives $(4,3)$ and $(5,3)$.

### Remaining profiles

- **$(4,4)$.** Each side has degree sequence $(2,1,1,1)$. A cycle is
  impossible because $C_4$ needs two degree-two vertices on each side. The
  three-component edge partitions are $(3,1,1)$ and $(2,2,1)$:
  $P_4+2K_2$, or two oppositely centered $P_3$'s plus $K_2$. The latter
  contains a crossing and is rejected. In $P_4+2K_2$, delete the central edge
  joining the unique degree-two row and column. The remaining four edges are
  the **only** all-distinct four-subset. All other deletions give repeated
  profiles and no crossing. Hence precisely that one determinant decides.
- **$(4,5)$.** Column degrees are all one and row degrees are $(2,1,1,1)$,
  forcing one row-centered $P_3+3K_2$. The two all-distinct four-subsets retain
  one or the other edge at the repeated row. All other deletions leave a
  repeated row and cannot be a size-four circuit. Thus exactly two determinants
  decide. The $(5,4)$ case is symmetric.
- **$(5,5)$.** Every degree is one, so the graph is a matching $5K_2$. Every
  deletion is an all-distinct four-matching; all five determinants must be
  nonzero.

This proves the template table exhaustively, independent of the census.

## 4. The tight $(1,2)$ family is a circuit

Choose a row triple

$$R=\{u_1,u_2,w\},$$

one active $C_1$ column $j_0$, and two $B$-circuit triples
$Z_1,Z_2$ containing $j_0$. Put the two active $C_2$ rows at $u_1,u_2$.
Let $c_R$ be the everywhere-nonzero relation on $R$ and let
$\delta_i$ be the relation on $Z_i$. Scale the row-$u_i$ relation by

$$\beta_i=-\frac{c_R(u_i)}{\delta_i(j_0)}.$$

Then the two overlap cells $(u_i,j_0)$ cancel, while every other coefficient
is nonzero. The support is exactly

$$S=\{(w,j_0)\}
\cup\bigl(\{u_1\}\times(Z_1\setminus\{j_0\})\bigr)
\cup\bigl(\{u_2\}\times(Z_2\setminus\{j_0\})\bigr).\tag{6}$$

Each summand is an axis relation, so their sum is a product relation. The graph
type depends only on $|Z_1\cap Z_2|$:

- 3: $C_4+K_2$, profile $(3,3)$;
- 2: $P_5+K_2$, profile $(3,4)$;
- 1 (only $j_0$): two row-centered $P_3$'s $+K_2$, profile $(3,5)$.

Section 3 proves all three templates have no fiber triple, no crossing
four-subset, and no all-distinct four-subset. Therefore every support (6) is a
genuine circuit. This also proves directly that a candidate cannot hide the
size-four crossing that would destroy minimality.

The unique $K_2$ component recovers $(w,j_0)$. The remaining component(s)
recover the two active rows and each incident $Z_i\setminus\{j_0\}$. Thus the
parameterization is injective (with the actual rows distinguishing the ordered
assignment when $Z_1\ne Z_2$).

## 5. Closed counts and their structural coincidence

Choose $R$ and $w$ in $3\binom n3$ ways.

### Shared triple: profile $(3,3)$

Choose $Z$ and $j_0\in Z$ in $3\binom n3$ ways. Hence

$$N_{33}=9\binom n3^2.$$

The transpose parameterization yields the same supports, not a second copy:
$C_4+K_2$ has a unique isolated edge and unique complementary rectangle.

This count coincides numerically with the frozen size-four crossing count for
a structural reason, not an accident. Both choose precisely
$(R,Z,w,j_0)$. The size-four crossing uses the two cross arms

$$\{w\}\times(Z\setminus\{j_0\})
\quad\cup\quad
(R\setminus\{w\})\times\{j_0\},$$

whereas the size-five family uses the complementary $2\times2$ rectangle
$(R\setminus\{w\})\times(Z\setminus\{j_0\})$ plus the opposite isolated edge
$(w,j_0)$.

### Intersection two: profile $(3,4)$

Fix $j_0$. Choose the other shared column ($n-1$ choices), an unordered pair
of distinct unique columns from the remaining $n-2$ points, and assign those
two unique columns to the two actual rows (2 choices). For fixed
$(R,w,j_0)$ this is

$$2(n-1)\binom{n-2}{2}=6\binom{n-1}{3}$$

choices. Thus

$$N_{34}=18n\binom n3\binom{n-1}3
=72\binom n3\binom n4.$$

### Intersection one: profile $(3,5)$

Fix $j_0$. Choose four other columns, partition them into two unordered pairs
(3 partitions), and assign the pairs to the two rows (2 choices). This gives
$6\binom{n-1}{4}$ choices per $(R,w,j_0)$, hence

$$N_{35}=18n\binom n3\binom{n-1}4.$$

Transposition gives $N_{43}=N_{34}$ and $N_{53}=N_{35}$. The derivation uses
only distinctness, no field values, so all five distinct profile counts are
field-independent.

## 6. Profile $(4,4)$: one determinant and a closed reduction

Let $S=P_4+2K_2$, with central edge $(u_1,v_1)$. Removing it gives the unique
all-distinct matching $Q$. By §3,

$$S\text{ is a circuit}\iff \det H_Q\ne0.\tag{7}$$

Conversely start with any independent four-matching $Q$. Choose an ordered pair
of distinct matching edges, say $(w_1,v_1)$ first and $(u_1,x_1)$ second, and
add $(u_1,v_1)$. The result is $P_4+2K_2$. Its unique central edge recovers
$Q$ and the ordered pair, so this parameterization is injective. There are
$4\cdot3=12$ ordered pairs. Since exactly
$M_4-N_{\rm all}^{(4)}(X,Y)$ all-distinct four-matchings are independent,

$$N_{44}^{(5)}=12\bigl(M_4-N_{\rm all}^{(4)}(X,Y)\bigr).$$

### Canonical $(2,2)$ decomposition

Write the four noncentral edges as

$$(w_1,v_1),\quad(w_2,v_2),\quad(u_1,x_1),\quad(u_2,x_2).$$

Use $C_1$ triples

$$R_1=\{u_1,u_2,w_1\}\text{ at }v_1,
\qquad R_2=\{u_1,u_2,w_2\}\text{ at }v_2,$$

and $C_2$ triples

$$Z_1=\{v_1,v_2,x_1\}\text{ at }u_1,
\qquad Z_2=\{v_1,v_2,x_2\}\text{ at }u_2.$$

Their four overlap cells form a square. Choose any nonzero scale for $R_2$;
then successively scale $Z_1,Z_2,R_1$ to cancel
$(u_1,v_2),(u_2,v_2),(u_2,v_1)$. Each division is by a nonzero triple-relation
coefficient. If the fourth overlap $(u_1,v_1)$ also canceled, the remaining
matching $Q$ would carry a nonzero relation, contradicting (7). It therefore
survives and the support is exactly $S$. Thus every accepted $(4,4)$ support
has

$$(k_1,k_2,T_1,T_2,r,c)=(2,2,6,6,4,3),$$

which demonstrates why the corrected law (1) is essential.

## 7. Profiles $(4,5),(5,4),(5,5)$

### The exact two-minor sum for $(4,5)$

Choose four rows $U$, five columns $V$, a repeated row $u\in U$, an unordered
pair $\{v_1,v_2\}\subset V$, and a bijection

$$f:U\setminus\{u\}\longrightarrow V\setminus\{v_1,v_2\}.$$

Let

$$Q_i=\{(u,v_i)\}\cup\operatorname{graph}(f),\qquad i=1,2.$$

Then the exact count is

$$N_{45}^{(5)}=
\sum_{U,V,u,\{v_1,v_2\},f}
\mathbf1[\det H_{Q_1}\ne0]\mathbf1[\det H_{Q_2}\ne0].\tag{8}$$

The graph data are recoverable from the support, so there is no quotient or
overcount. The number of unfiltered candidates is

$$240\binom{|X|}{4}\binom{|Y|}{5}$$

because $4\binom52 3!=240$ for fixed $U,V$. The $(5,4)$ formula is the
transpose of (8).

### Canonical $(2,2)$ decomposition for $(4,5)$

Write the repeated-row edges as $(u_1,x_1),(u_1,x_2)$ and the three isolated
edges as $(u_2,x_3),(w_1,v_1),(w_2,v_2)$. Use $C_1$ triples

$$\{u_1,u_2,w_1\}\text{ at }v_1,
\qquad\{u_1,u_2,w_2\}\text{ at }v_2,$$

and a $C_2$ triple $\{v_1,v_2,x_3\}$ at $u_2$. Scale these three relations to
cancel $(u_2,v_1),(u_2,v_2)$. The required coefficients of the row-$u_1$
relation at $v_1,v_2$ are then fixed as their negatives. Since
$b_{x_1},b_{x_2}$ are independent, there is a unique pair of coefficients at
$x_1,x_2$ that completes a $B$ relation; these cancel the remaining two
overlap cells. The resulting relation is nonzero and supported inside $S$.
Because (8) makes $S$ a circuit, no intended survivor coefficient can vanish,
or a proper subset would carry that relation. Thus the row-$u_1$ relation has
support four, and every accepted support has

$$(k_1,k_2,T_1,T_2,r,c)=(2,2,6,7,4,4).$$

The transpose has $(T_1,T_2)=(7,6)$.

### The exact five-minor sum for $(5,5)$

Choose five rows $U$, five columns $V$, and a bijection $f:U\to V$. Then

$$N_{55}^{(5)}=
\sum_{U,V,f}\prod_{u\in U}
\mathbf1\left[
\det H_{\operatorname{graph}(f|_{U\setminus\{u\}})}\ne0
\right].\tag{9}$$

Equivalently the five product columns have rank four and their unique relation
has full support. The weaker factor-rank inequality

$$\rho_A+\rho_B\le5$$

is necessary but not sufficient: here every candidate has
$\rho_A=\rho_B=2$, so all candidates satisfy $4\le5$, including those with a
zero cofactor. At $n=5,p=13$, all 120 matchings satisfy the factor-rank
inequality but only 60 satisfy (9).

### Canonical $(3,2)$ decomposition for $(5,5)$

Let the matching be $(u_i,v_i)$, $i=1,\dots,5$, with full-support relation

$$\sum_{i=1}^5\gamma_i a_{u_i}\otimes b_{v_i}=0,
\qquad\gamma_i\ne0.$$

Choose $a_{u_1},a_{u_2}$ as a basis and let $\ell_1,\ell_2$ be its dual
coordinate functionals. Applying $\ell_r\otimes\mathrm{id}$ to the relation
gives two $B$ relations

$$\sum_i\ell_r(a_{u_i})\gamma_i b_{v_i}=0,
\qquad r=1,2.$$

Use them as the two nonzero rows of $C_2$. We have
$\ell_1(a_{u_2})=0$ and $\ell_2(a_{u_1})=0$. Every other coordinate is
nonzero: if, for example, $\ell_1(a_{u_i})=0$ for $i\ne2$, then $a_{u_i}$ is
parallel to $a_{u_2}$, contradicting spark three. Thus each $C_2$ row has
support four, for $T_2=8$.

Set $C_1=\Gamma-C_2$. Its columns at $v_1,v_2$ vanish. For $i=3,4,5$, the
column at $v_i$ is the three-term relation

$$\gamma_i\bigl(a_{u_i}-\ell_1(a_{u_i})a_{u_1}
-\ell_2(a_{u_i})a_{u_2}\bigr)=0,$$

with all three coefficients nonzero. Hence $C_1$ has three minimum-support
columns, $T_1=9$. The six cells
$\{u_1,u_2\}\times\{v_3,v_4,v_5\}$ cancel, while the five matching cells
survive. Therefore

$$(k_1,k_2,T_1,T_2,r,c)=(3,2,9,8,6,6),$$

and (1) gives $9+8-6-6=5$. Transposition gives a $(2,3)$ construction.

## 8. Exact census and set-equality evidence

The standalone instrument used exact Python integers, actual columns
$(1,y,x,xy)$, exact Gauss–Jordan ranks/determinants, and no frozen import. For
each of six configurations it swept every five-subset by two routes:

1. actual rank four plus all five deletion ranks four;
2. graph template, crossing exclusion, and only the prescribed all-distinct
   minors.

The full support sets were equal in every row. The independently constructed
$(1,2)/(2,1)$, P44, P45/P54, and P55 support sets also equaled their measured
profiles exactly. Every surviving P44/P45/P54/P55 support was reconstructed by
its canonical decomposition.

### $n=4$: all $\binom{16}{5}=4368$ supports

| prime | total | $(3,3)$ | $(3,4)$ | $(4,3)$ | $(4,4)$ | size-4 $N_{\rm all}$ |
|---:|---:|---:|---:|---:|---:|---:|
| 7 | 912 | 144 | 288 | 288 | 192 | 8 |
| 11 | 960 | 144 | 288 | 288 | 240 | 4 |
| 13 | 864 | 144 | 288 | 288 | 144 | 12 |

In every row, $N_{44}=12(24-N_{\rm all})$. All 432 oriented $(1,2)$
parameters had exact axis-kernel coefficients and were injective; the transpose
set agreed, including the common 144 profile-$(3,3)$ supports.

### $n=5$: all $\binom{25}{5}=53130$ supports

| $p$ | total | 33 | 34 | 35 | 43 | 44 | 45 | 53 | 54 | 55 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 7 | 17124 | 900 | 3600 | 900 | 3600 | 5760 | 720 | 900 | 720 | 24 |
| 11 | 17560 | 900 | 3600 | 900 | 3600 | 6000 | 800 | 900 | 800 | 60 |
| 13 | 17880 | 900 | 3600 | 900 | 3600 | 6192 | 864 | 900 | 864 | 60 |

The size-four all-distinct counts entering P44 were respectively 120, 100,
84, so $12(600-N_{\rm all})$ gives 5760, 6000, 6192 exactly. P45 filtered
1200 graph candidates to 720/800/864; P55 filtered 120 matchings to 24/60/60.
All 5400 oriented $(1,2)$ parameters per prime had exact verified relation
coefficients and multiplicity one.

These rows expose the exact boundary: profiles with a 3-coordinate are
field-free; P44/P45/P55 are determinant-incidence counts and vary with the
field.

## 9. Fail-loud controls

- **Known-true ACCEPTs.** One support from each of $(3,3),(3,4),(3,5)$ had
  verified $(1,2)$ coefficients, rank four, all five nonzero deletion minors,
  and acceptance by both routes.
- **Inserted crossing REJECT.** A known four-edge crossing plus a fifth cell
  retained its zero determinant as a proper subset; both direct and structural
  routes rejected it.
- **Factor-rank-only REJECT.** At $n=5,p=13$, the inequality accepted all 120
  matchings, while the five-minor predicate accepted exactly 60.
- **Non-tensor label REJECT.** A label-valid $(3,4)$ support remained accepted
  structurally after one actual product column was replaced by a duplicate of
  another; actual rank/proper-subset checks rejected the corrupted columns.
- **Duplicate factor REJECT.** Cloning one factor column dropped measured spark
  to two and exposed exactly the planted dependent pair $(0,5)$.
- **Kronecker/concat guard.** Both true $2\times2$ Kronecker coordinates and
  block concatenation have dimension four, so dimension alone cannot catch the
  bug. Exact coordinates $[3,4,6,8]$ rejected concatenation $[1,2,3,4]$.
- **Wrong support law REJECT.** A certified P44 decomposition gave
  $(6,6,4,3)$: (1) returned five and unrestricted `-2c` returned six.

All 170 aggregate assertions passed. CPU was 8.591/600 seconds, wall time
9.511/900 seconds, nice value 15, thread caps 1, exact arithmetic only, and
bytecode disabled.

## 10. Scope, defects, and verdict

**Certified:** the complete size-five support/template theorem for two 2-row
GRS factors; the corrected general support-cardinality identity; the exhaustive
irreducible tight-shape result; genuine-circuit and injectivity proofs for the
$(1,2)/(2,1)$ families; all field-free closed counts; P44's closed reduction;
P45/P54 and P55 exact determinant sums; and canonical decomposition shapes for
every accepted profile.

**Not claimed:** a field-free closed polynomial for P45/P55; size six or
larger; arbitrary higher-row factors; or spark-two factors. The simple
$\rho_A+\rho_B\le5$ inequality is explicitly not promoted as a sufficient
P55 test, and (2) is explicitly not promoted beyond its tight/full-cancel
hypotheses.

No run defect, failed assertion, prereg amendment, or parameter amendment
occurred. Every supplied anchor was reproduced, all three primes passed full
set equality for both $n=4,5$, every canonical decomposition was reconstructed,
and every plant fired. The exact supported verdict is **FROZEN-CERTIFIED**.
