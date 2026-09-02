# THEOREM H-3ROW-DETERMINANTAL
## Exact residual circuits for a 2-row by 3-row GRS product

Run `20260902T025427Z_1b1f202c_6247a0994bce`, agent `RsPe3dH2`,
2026-09-02 UTC. Preregistration:
`prereg/H_3ROW_DETERMINANTAL_PREREG_2026-09-02.md`, path-scoped commit
`a01f1e598a29ee81daac8bd7d811f8d068dd0ad9`, SHA-256
`dae23f1012f9509ed37d7e8b85e6852d44562e6b1e01977fd3a77f985067062c`,
with byte-identical in-run copy `pre_statement.md`. There was no preregistration
amendment and no parameter-dependent computation before the commit and run
initialization.

The verdict is **FROZEN-CERTIFIED**. The size-five profile-$(5,5)$ and size-six
profiles $(4,5),(5,5)$ left open by the predecessor are classified below by
Möbius incidence, with injective exact PGL sums. Every residual size-seven
channel has a necessary-and-sufficient deletion-minor predicate. The remaining
field-dependent size-seven sums are exact finite sums, not falsely promoted to
field-free polynomials.

## 1. Evaluation matrix and the exact circuit principle

Let $X,Y$ be finite sets of distinct field elements and put

\[
a_x=(1,x)^{\mathsf T},\qquad b_y=(1,y,y^2)^{\mathsf T},
\]

\[
h(x,y)=a_x\otimes b_y=(1,y,y^2,x,xy,xy^2)^{\mathsf T}.
\]

For a support $S$, let $E_S$ be the matrix with columns $h(x,y)$ for
$(x,y)\in S$. A indices are row vertices and B indices are column vertices.
Here $(d_A,d_B)=(3,4)$, so every nonfiber size-5--7 circuit has

\[
|\pi_A S|\ge3,\quad |\pi_BS|\ge4,
\quad \deg_A\le3,
\quad \deg_B\le2. \tag{1}
\]

The following is the exact general principle used in this run.

**Circuit principle.** For $m=|S|\le7$,

\[
S\text{ is a circuit}
\iff \operatorname{rank}E_S=m-1
\text{ and }\operatorname{rank}E_{S\setminus\{e\}}=m-1
\quad(e\in S). \tag{2}
\]

The forward implication is the definition of a linear-matroid circuit. For the
reverse implication, every deletion is independent; independence is hereditary,
so every proper subset is independent, while $S$ is dependent with relation
space of dimension one. In the present ambient dimension six this specializes
to:

- $m=5$: rank four and all four-deletions rank four;
- $m=6$: one dependency condition plus all five-deletions rank five;
- $m=7$: dependence is automatic and all seven six-deletion determinants must
  be nonzero.

Thus “minimal graph” alone is not a theorem when field-dependent minors remain;
the graph and the actual evaluation minors are both load-bearing.

## 2. The six-column determinant

For an ordered six-set $T=\{(x_i,y_i):1\le i\le6\}$ write

\[
D_6(T)=\det[h(x_i,y_i)]_{i=1}^6.
\]

Laplace expansion along the first three rows gives

\[
\boxed{D_6(T)=
\sum_{I\subset[6],\ |I|=3}(-1)^{\sum_{i\in I}i}
\left(\prod_{j\notin I}x_j\right)V(y_I)V(y_{I^c})}, \tag{3}
\]

where $V(y_i,y_j,y_k)$ is the ordered $3\times3$ Vandermonde determinant.
This is an explicit polynomial in the evaluation points, not a rank oracle.
The run compared (3) with a separately implemented exact $6\times6$
determinant on every six-subset: 8,008 per $n=4$ configuration and 177,100
per $n=5$ configuration, at all three primes, with exact equality.

There is also a useful left-kernel interpretation. The square matrix is
singular exactly when a nonzero polynomial

\[
P(x,y)=f(y)+xg(y),\qquad \deg f,\deg g\le2, \tag{4}
\]

vanishes at all six cells. Sorting by the number $q$ of distinct B values gives
four reductions.

### 2.1 Three B values: degrees $(2,2,2)$

Each B column is paired with two distinct A columns. The three B columns form
a basis of the B space and each A pair forms a basis of the A space. The six
product columns are therefore a tensor basis, so $D_6\ne0$.

### 2.2 Four B values: degrees $(2,2,1,1)$

At either repeated B value, the two distinct A evaluations force both $f$ and
$g$ to vanish. Hence $f,g$ are multiples of the quadratic having the two
repeated B values as roots. At the two singleton B values, (4) then has a
nonzero solution exactly when the singleton edges have the same A endpoint:

\[
\boxed{D_6=0\iff x_{\rm singleton,1}=x_{\rm singleton,2}}. \tag{5}
\]

### 2.3 Five B values: degrees $(2,1,1,1,1)$

Let $y_0$ be repeated at $x_a\ne x_b$, and let $R$ be the four singleton-B
edges. Subtracting the two equations at $y_0$ gives $g(y_0)=0$, and then
$f(y_0)=0$. Thus both have the factor $y-y_0$, and their linear quotients
annihilate $R$. If the repeated columns occupy zero-indexed ordered positions
$a<b$, direct expansion gives the exact identity

\[
\boxed{D_6(S)=(-1)^{b-a-1}(x_b-x_a)
\prod_{(x,y)\in R}(y-y_0)D_4(R)}, \tag{6}
\]

where

\[
D_4(R)=\det[1,y,x,xy]_{(x,y)\in R}. \tag{7}
\]

All prefactors in (6) are nonzero, so the single dependency condition is
$D_4(R)=0$. The run checked (6) on all 31,250 applicable six-sets at each of
$p=7,11,13$; the observed sign was exactly $(-1)^{b-a-1}$ in all 15 possible
position classes.

### 2.4 Six B values

Formula (3) is the irreducible predicate needed here. Equivalently,
$D_6=0$ iff the six paired points lie on a possibly degenerate rational graph

\[
x=-f(y)/g(y),\qquad \deg f,\deg g\le2, \tag{8}
\]

with the polynomial equation (4) interpreted at poles. No uniqueness or
nondegeneracy of (8) is asserted without additional hypotheses.

## 3. Complete size-five theorem

A size-five circuit is exactly one of two types:

1. a profile-$(3,4)$ crossing;
2. a profile-$(5,5)$ matching lying on a unique nondegenerate Möbius graph.

The first type is the certified Theorem-X family and has count

\[
12\binom{|X|}{3}\binom{|Y|}{4}. \tag{9}
\]

It remains to prove the second type and exclude all intermediate profiles.

### 3.1 Five distinct B coordinates force a Möbius graph

Let $y_1,\ldots,y_5$ be distinct and define barycentric weights

\[
w_j=\left(\prod_{k\ne j}(y_j-y_k)\right)^{-1}.
\]

The standard partial-fraction/Vandermonde identities give

\[
\sum_jw_jy_j^t=0\quad(0\le t\le3),
\]

and therefore

\[
\ker[1;y;y^2]=\operatorname{span}\{(w_j),(w_jy_j)\}. \tag{10}
\]

A product-column relation $\lambda$ is equivalent to

\[
[1;y;y^2]\lambda=0,
\qquad [1;y;y^2]D_x\lambda=0. \tag{11}
\]

By (10), there are $\alpha,\beta,\gamma,\delta$ such that

\[
\lambda_j=w_j(\alpha+\beta y_j),
\qquad
x_j(\alpha+\beta y_j)=\gamma+\delta y_j. \tag{12}
\]

For a circuit, $\lambda$ has full support. Hence the first linear form in
(12) is nonzero at every selected $y_j$ and

\[
x_j=M(y_j):=\frac{\gamma+\delta y_j}{\alpha+\beta y_j}. \tag{13}
\]

The map is nondegenerate: a degenerate ratio is constant, contradicting five
distinct A endpoints. Conversely, (13) makes
$\lambda_j=w_j(\alpha+\beta y_j)$ a full-support product relation. Every four
matching edges are independent: the B Vandermonde kernel is then
one-dimensional, and (11) would force all four $x_j$ equal. Thus the five-set
is a circuit. Three paired points determine $M$, proving uniqueness.

For $M\in\mathrm{PGL}(2,F)$ define

\[
K_M=\{y\in Y:M(y)\in X\},\qquad k_M=|K_M|.
\]

A five-subset of $K_M$ gives one matching, and uniqueness prevents double
counting. Therefore

\[
\boxed{N_{55}^{(5)}(X,Y)=
\sum_{M\in\mathrm{PGL}(2,F)}\binom{k_M}{5}}. \tag{14}
\]

### 3.2 Four distinct B coordinates force a crossing

Now one B value occurs twice and three are singleton values. The kernel of the
B Vandermonde matrix is generated by the duplicate-column difference and the
unique full relation on the four distinct B columns. In a full product
relation, the three singleton coefficients are nonzero scalar multiples of
the latter relation. Applying $D_x$ in (11) shows their three A values are all
equal. The two cells at the repeated B value have distinct A values, and the
degree bound prevents either from being the singleton center. The support is
exactly a profile-$(3,4)$ crossing. This excludes every other size-five
profile analytically.

## 4. Complete size-six residual theorem

### 4.1 The unified field-free $q=4$ family

By (5), the two singleton-B edges must share an A row $r_0$. The six-set is
minimal exactly when $r_0$ is not an endpoint of either doubled-B pair.
Otherwise the two singleton edges, the incident doubled-B edge, and the other
doubled-B pair form a proper size-five crossing. If $r_0$ avoids both pairs,
there is no crossing and no profile-$(5,5)$ matching proper circuit, so every
five-deletion is independent.

For profile $(p,4)$, $p=3,4,5$, let

\[
h_q=\sum_{t=0}^q(-1)^t\binom qt\binom{q-t}{2}^2. \tag{15}
\]

After choosing the unique row $r_0$, the two labelled doubled B vertices choose
two pairs among the other $p-1$ A vertices whose union must cover all of them;
(15) is inclusion-exclusion for that choice. The injective count is

\[
\boxed{N_{p4}^{(6)}=6p\binom{|X|}{p}\binom{|Y|}{4}h_{p-1}}. \tag{16}
\]

Here $h_2=1,h_3=6,h_4=6$, recovering the coefficients 18, 144, and 180 of the
three predecessor families.

### 4.2 Profile $(4,5)$: $P_5+2K_2$

Let $R$ be the four singleton-B edges. Equation (6) makes $D_4(R)=0$ the one
condition carrying dependence. Circuit minimality rules out a degenerate
bilinear graph, since its three constant A values together with the repeated-B
pair contain a crossing. Thus $R$ is the graph of a unique nondegenerate
$M$, uses four A values, and the repeated B vertex joins two of them. The graph
is $P_5+2K_2$.

No deletion can be a size-five matching, and direct inspection of the deletion
graphs shows none is a crossing. By §3 all five-deletions are therefore
independent. Hence

\[
\boxed{S\text{ is a circuit}\iff D_4(R)=0}. \tag{17}
\]

Choose $M$, four values in $K_M$, a new B value outside those four, and two of
the four graph edges to receive the common new B neighbor. The graph recovers
all these parameters, so

\[
\boxed{N_{45}^{(6)}(X,Y)=
6(|Y|-4)\sum_M\binom{k_M}{4}}. \tag{18}
\]

### 4.3 Profile $(5,5)$: $P_4+3K_2$

Again let $R$ be the four singleton-B edges. Minimal dependence forces $R$ to
be a nondegenerate Möbius matching. Exactly one endpoint of the repeated-B
pair belongs to the four A values of $R$, and the other is new. Thus the graph
is $P_4+3K_2$; its central edge joins the unique degree-two A and B vertices.

Deleting the central edge leaves a five-edge matching $Q$. All other deletions
have a profile excluded by §3 and are automatically independent. Consequently

\[
\boxed{S\text{ is a circuit}
\iff D_4(R)=0\text{ and }Q\text{ is not a Möbius five-circuit}}. \tag{19}
\]

If $R$ is on $M$, the second condition is simply that the new A endpoint at the
repeated B value $y_0$ is not $M(y_0)$. Choose $M$, four values in $K_M$, the
distinguished edge of $R$ incident to the central edge (4 choices), a new B
value, and a new A value outside $M(R)$. For each fixed four-set, there are
$(|X|-4)(|Y|-4)$ new A/B pairs, of which $k_M-4$ extend the same Möbius graph
and must be excluded. The parameters are graph-recoverable, giving

\[
\boxed{N_{55}^{(6)}(X,Y)=
4\sum_M\binom{k_M}{4}
\bigl((|X|-4)(|Y|-4)-(k_M-4)\bigr)}. \tag{20}
\]

The omitted condition in (19) is genuinely necessary. On $X=Y=\{1,\ldots,5\}$,
take four identity-graph edges, add a new B value to one old A endpoint, and
also add its identity A endpoint. Then $D_4(R)=0$, but deleting the central
edge leaves the five-edge identity Möbius circuit. The wrong predicate
`D4(R)=0` overaccepted 120, 200, and 120 supports at $p=7,11,13$ respectively;
(19) had exact set equality.

## 5. Complete size-seven deletion-minor theorem

A seven-set is dependent solely by the ambient dimension. There is no extra
“single vanishing determinant”: the one condition carrying dependence is
$7>6$. Circuitness is exactly the nonvanishing of all six-deletion
determinants. Equations (3), (5), and (6) reduce this to the following explicit
predicate according to the B-degree sequence.

### 5.1 $q=4$: B degrees $(2,2,2,1)$

Let $(r_0,y_0)$ be the unique singleton-B edge. Deleting it leaves three B
basis vectors, each paired with an A basis, and is automatically independent.
Deleting an edge at a doubled B vertex leaves B degrees $(2,2,1,1)$; by (5),
that minor vanishes exactly when the remaining endpoint of the affected pair
is $r_0$. Hence

\[
\boxed{S\text{ is a circuit}\iff
r_0\text{ occurs in none of the three doubled-B pairs}}. \tag{21}
\]

This channel is completely field-free. For profile $(p,4)$, $3\le p\le7$,
let

\[
g_q=\sum_{t=0}^q(-1)^t\binom qt\binom{q-t}{2}^3. \tag{22}
\]

The three labelled doubled B vertices choose three pairs among the $p-1$ A
vertices other than $r_0$, with union all $p-1$. Thus

\[
\boxed{N_{p4}^{(7)}=4p\binom{|X|}{p}\binom{|Y|}{4}g_{p-1}}. \tag{23}
\]

For the swept profiles, $g_2=1,g_3=24,g_4=114$, giving coefficients 12, 384,
and 2,280 per selected A/B vertex sets.

### 5.2 $q=5$: B degrees $(2,2,1,1,1)$

Let the three singleton-B edges have A endpoints $r_1,r_2,r_3$. Deleting one
of them leaves a $q=4$ six-set whose two singleton endpoints are the other two
$r$'s. By (5), all three such minors are nonzero exactly when

\[
r_1,r_2,r_3\text{ are pairwise distinct}. \tag{24}
\]

Deleting an edge from a doubled-B pair leaves one repeated B value and four
singleton edges: the retained endpoint of the affected pair together with the
original three singleton edges. By (6), its determinant is nonzero exactly
when the corresponding $D_4$ is nonzero. There are two endpoints in each of
two pairs, hence four conditions:

\[
D_4(\{\text{three singleton edges}\}\cup\{e\})\ne0
\quad\text{for every endpoint edge }e\text{ of a doubled B vertex}. \tag{25}
\]

Conditions (24)--(25) are necessary and sufficient. They apply verbatim to all
component signatures and profiles with $q=5$; field dependence lives precisely
in the four bilinear minors.

### 5.3 $q=6$: B degrees $(2,1,1,1,1,1)$

Deleting any singleton-B edge leaves the repeated pair and four singleton
edges, so require the five corresponding $D_4$ values from (6) to be nonzero.
Deleting either repeated-B edge leaves six distinct B values, so require the
two determinants (3) to be nonzero. These seven conditions are necessary and
sufficient.

### 5.4 $q=7$: B degrees $(1,1,1,1,1,1,1)$

Every deletion has six distinct B values. Require all seven determinants (3)
to be nonzero. This is again necessary and sufficient.

This proves the general size-seven predicate for every possible B projection
$q=4,5,6,7$, not only the $n\le5$ profiles swept in-run.

## 6. Exact finite corroboration

`run_determinantal.py` independently computed direct circuit sets from actual
rank/deletion tests and predicted sets from §§2--5. It asserted full support-set
equality overall, by profile, and for every measured size-seven
component/degree signature at each of six configurations:

\[
n\in\{4,5\},\qquad p\in\{7,11,13\}.
\]

### 6.1 Census totals

For $n=4$, all three primes gave identical sets:

| size | universe | profiles | circuits |
|---:|---:|---|---:|
| 5 | 4,368 | $(3,4):48$ | 48 |
| 6 | 8,008 | $(3,4):72,(4,4):144$ | 216 |
| 7 | 11,440 | $(3,4):48,(4,4):384$ | 432 |

For $n=5$:

| size | GF(7) | GF(11) | GF(13) |
|---:|---|---|---|
| 5 | $(3,4):600,(5,5):6$ | $(3,4):600,(5,5):10$ | $(3,4):600,(5,5):6$ |
| 6 | $(3,4):900,(4,4):3600,(4,5):720,(5,4):900,(5,5):360$ | $(3,4):900,(4,4):3600,(4,5):600,(5,4):900,(5,5):200$ | $(3,4):900,(4,4):3600,(4,5):504,(5,4):900,(5,5):216$ |
| 7 | $(3,4):600,(3,5):5400,(4,4):9600,(4,5):23760,(5,4):11400,(5,5):12600$ | $(3,4):600,(3,5):5400,(4,4):9600,(4,5):25200,(5,4):11400,(5,5):15000$ | $(3,4):600,(3,5):5400,(4,4):9600,(4,5):26352,(5,4):11400,(5,5):15768$ |

Thus totals at sizes 5,6,7 were respectively
$606/6480/63360$, $610/6200/67200$, and $606/6120/69120$ across the three
primes.

At size seven there were 17 label-invariant component/degree signatures:
$1,2,2,4,3,5$ signatures in profiles
$(3,4),(3,5),(4,4),(4,5),(5,4),(5,5)$. The reduced predicate equalled the
direct set separately on every one. Profiles with $q=4$ had the exact
field-free populations from (23): 600, 9,600, and 11,400 at $n=5$. The
$q=5$ profile $(3,5)$ also happened to be set-stable across the three primes,
but no field-free formula is inferred from that observation; it is governed
by (24)--(25).

### 6.2 PGL incidence spectra and formula checks

For $X=Y=\{1,\ldots,5\}$, the exact $k_M$ histograms were

| prime | $k=0$ | 1 | 2 | 3 | 4 | 5 | $|\mathrm{PGL}(2,p)|$ |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 7 | 0 | 0 | 60 | 180 | 90 | 6 | 336 |
| 11 | 60 | 200 | 700 | 300 | 50 | 10 | 1,320 |
| 13 | 162 | 594 | 1,044 | 324 | 54 | 6 | 2,184 |

The independently constructed support sets from (14), (18), and (20) were
injective and equalled the direct sets exactly. Their counts were:

| prime | size-5 $(5,5)$ | size-6 $(4,5)$ | size-6 $(5,5)$ |
|---:|---:|---:|---:|
| 7 | 6 | 720 | 360 |
| 11 | 10 | 600 | 200 |
| 13 | 6 | 504 | 216 |

The support sets genuinely vary. For example, the size-five $(5,5)$ pairwise
intersections are all 2; the size-six $(4,5)$ intersections are all 216; and
the size-six $(5,5)$ intersections are 56 for $7/11$, 72 for $7/13$, and 56
for $11/13$. These changes are exactly the PGL-incidence mechanism, not census
noise.

For size-seven profiles $(4,5)$ and $(5,5)$ the pairwise intersections were
respectively $19{,}512/20{,}664/21{,}744$ and
$9{,}360/10{,}160/11{,}928$ for prime pairs $7/11$, $7/13$, $11/13$.
This confirms the field dependence already explicit in (25).

## 7. Controls, defects, and runtime

All 261 aggregate assertions passed. Positive plants were an actual crossing,
the five-edge identity Möbius matching, and a genuine seven-circuit with every
six-deletion independent. Reject plants were a larger set containing the
identity Möbius circuit, a nonminimal ambient-dependent seven-set, the
profile-$(5,5)$ wrong-minor example from §4.3, and a label-valid
profile-$(4,5)$ circuit after one actual product-column coordinate was
corrupted. A true Kronecker column had six coordinates while the planted block
concatenation had five, so the concat impostor was rejected.

Two defects were preserved and fully disclosed:

1. `broken_results.json`: the first invocation stopped before any census at a
   combined startup bytecode guard, even though its own detail recorded both
   operands true. Separate standalone and in-module probes passed. The guard
   was split into two named assertions; the failed run consumed 0.0001 CPU
   seconds and performed no parameter computation.
2. `broken_results_1788318130.json`: the second invocation incorrectly demanded
   a positive $q=5$ factorization population at $n=4$, where only four B values
   exist. The fix requires zero instances at $n<5$ and positive exhaustive
   coverage at $n\ge5$. The full six-configuration matrix was rerun.

The corresponding records are `defect_001_startup_guard.json` and
`defect_002_n4_factor_scope.json`. Neither changed the preregistered theorem or
parameters, so there was no preregistration amendment.

The successful all-family rerun used exact Python integers modulo $p$, stdlib
only, true row-major Kronecker columns, `sys.dont_write_bytecode=True`,
`PYTHONDONTWRITEBYTECODE=1`, `nice -n 15`, and all four declared thread counts
pinned to one. It consumed 83.035 CPU seconds and 89.815 wall seconds, below the
hard 900/1200-second caps.

## 8. Certified scope and remaining work

**Certified:** the complete size-five classification and PGL sum (14); the
six-column reductions (3), (5), (6); the exact profile-$(4,5)$ and $(5,5)$
size-six predicates and injective sums (18), (20); the unified field-free
$q=4$ size-six formula (16); and the complete size-seven deletion-minor
predicate (21), (24)--(25), and its $q=6,7$ extensions, including the
field-free $q=4$ formula (23).

**Not claimed:** a simpler closed evaluation of the general PGL incidence sums
for arbitrary $X,Y$; a field-free polynomial for the size-seven $q=5,6,7$
minor sums; uniqueness of the rational quadratic representation (8); circuit
classification above size seven or for higher factor row ranks. Those are
further simplification/generalization questions, not gaps in this gate's
explicit circuit predicates.
