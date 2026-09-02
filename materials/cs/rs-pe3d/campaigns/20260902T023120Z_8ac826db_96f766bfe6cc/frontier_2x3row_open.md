# Exact opening of the 2-row by 3-row GRS circuit frontier
## gate `H-2ROW-COMPLETE-3ROW-OPEN`, run `20260902T023120Z_8ac826db_96f766bfe6cc`

This is the higher-row companion to `theorem_2row_complete.md`. Its status is
**FROZEN-INCONCLUSIVE**: the census is exhaustive on the preregistered finite
instances, and several channels below are proved field-free, but sizes 5--7
are not analytically complete for arbitrary point sets. No census count is
promoted to a general formula without a proof.

Preregistration commit:
`4486a168f8cc11a07f3e1585531e768502ca3408`; preregistration SHA-256:
`18f60f53b292e873f31d01231f8634d24250c37e237eabb76c2859c402ee60b1`.
The in-run `pre_statement.md` is byte-identical and there was no amendment.

## 1. Setting and universal constraints

Put

\[
a_i=(1,x_i)^{\mathsf T},\qquad
b_j=(1,y_j,y_j^2)^{\mathsf T},\qquad
h_{ij}=a_i\otimes b_j=(1,y_j,y_j^2,x_i,x_i y_j,x_i y_j^2)^{\mathsf T}.
\]

The evaluation points are distinct within each factor. Thus
$(d_A,d_B)=(3,4)$, product ambient dimension is six, and every circuit has size
at most seven by the rank bound. The product spark is three, so possible sizes
are exactly in the interval 3--7.

A nonfiber circuit has projection floor

\[
|\pi_A S|\ge3,\qquad |\pi_BS|\ge4. \tag{1}
\]

To see this, if fewer than three A columns occur, those A columns are
independent. Dual functionals isolate each A coordinate in a relation and
force each active row section to be a B relation, hence to contain at least
four cells. Circuit minimality then forces a single row, which is a B fiber.
The argument with the factors reversed proves the other inequality.

For a circuit of size 5--7, every factor fiber is a proper subset. Therefore,
with A indices drawn as row vertices and B indices as column vertices,

\[
\deg_A\le3=d_B-1,\qquad \deg_B\le2=d_A-1. \tag{2}
\]

The side labels in (2) matter: the opposite factor supplies the forbidden
fiber.

The already-closed layers are:

- size 3: exactly A-factor triples at one B index, profile $(3,1)$, count
  $n\binom n3$;
- size 4: exactly B-factor quadruples at one A index, profile $(1,4)$, count
  $n\binom n4$.

There is no nonfiber size-four channel. A crossing would have size
$d_A+d_B-2=5$, and the certified near-minimum necessary condition excludes an
all-distinct size-four circuit because $d_A+d_B=7>d+3=6$.

## 2. Exact circuit predicates used throughout

For a support $S$, let $H_S$ be the six-by-$|S|$ matrix of its actual product
columns. These predicates are definitions specialized to the ambient rank and
were evaluated by exact modular Gauss--Jordan elimination, never by labels:

- $|S|=5$: $S$ is a circuit iff $\operatorname{rank}H_S=4$ and every
  four-deletion has rank four;
- $|S|=6$: $S$ is a circuit iff $\operatorname{rank}H_S=5$ and every
  five-deletion has rank five;
- $|S|=7$: $S$ is automatically dependent, but is a circuit iff every
  six-deletion has rank six.

Thus the unresolved channels below are still given exact finite determinant
predicates. “OPEN” means that the corresponding determinant sum has not been
converted into a field-free parameterization and closed count, not that the
finite predicate is ambiguous.

## 3. Size five

The complete Theorem-X crossing converse applies at
$d_A+d_B-2=5$. Every profile-$(3,4)$ five-circuit is therefore a crossing

\[
S(T,i_0,Z,j_0)=((T\setminus\{i_0\})\times\{j_0\})
\cup(\{i_0\}\times(Z\setminus\{j_0\})), \tag{3}
\]

where $|T|=3$, $|Z|=4$, and conversely every such support is a circuit. Its
parameters are recovered from the unique degree-three A vertex and unique
degree-two B vertex, so the parameterization is injective and

\[
\boxed{N_{34}^{(5)}=12\binom n3\binom n4}. \tag{4}
\]

The only additional measured profile was $(5,5)$, necessarily a five-edge
matching. Its exact predicate is

\[
\operatorname{rank}H_S=4
\quad\text{and}\quad
\operatorname{rank}H_{S\setminus e}=4\quad(e\in S). \tag{5}
\]

Equivalently all $5\times5$ minors of $H_S$ vanish while every four-deletion
has some nonzero $4\times4$ minor. This is a determinantal channel, not a
field-free graph channel. At $n=5$ it has six supports at each tested prime,
but only two supports are common: each prime has four supports absent at the
other. The equality of the two totals is therefore not a general law.

## 4. Three proved field-free size-six families

For compactness, $A_U=\operatorname{span}\{a_i:i\in U\}$ and similarly for
$B_V$. Any two distinct A columns form a basis of the two-dimensional A space;
any three distinct B columns form a basis of the three-dimensional B space.
When $V,W$ are disjoint B pairs,

\[
\dim(B_V\cap B_W)=1, \tag{6}
\]

and the vector on that line has both coordinates nonzero in either pair basis.
Otherwise three distinct B columns would be dependent. Likewise, a third A
column has two nonzero coordinates in any basis of two distinct A columns.
These elementary MDS facts prove the following families over every field with
the required distinct points.

### 4.1 Profile $(3,4)$: $C_4$ plus an A-centered $P_3$

Choose three A vertices and four B vertices. Choose two A and two B vertices
for a complete $2\times2$ rectangle; join the remaining A vertex to the two
remaining B vertices. The rectangle spans
$A\otimes B_V$ (dimension four), the $P_3$ spans $a_w\otimes B_W$
(dimension two), and their intersection is
$a_w\otimes(B_V\cap B_W)$ (dimension one). Hence the six columns have rank
five.

If a rectangle edge is deleted, the common-line vector has all four nonzero
coordinates in the rectangle tensor basis, so it is outside the resulting
three-dimensional coordinate hyperplane. If a $P_3$ edge is deleted, its
remaining B column is outside $B_V$. Every five-deletion therefore has rank
five. The support is a circuit.

The two connected components distinguish the rectangle and $P_3$, so the
parameters are recovered from the support. The exact count is

\[
\boxed{N_{34}^{(6),\mathrm{ff}}=18\binom n3\binom n4}. \tag{7}
\]

### 4.2 Profile $(4,4)$: $P_5$ plus an A-centered $P_3$

Choose four A and four B vertices. The $P_5$ uses three A vertices and two B
vertices, with its degree-two A center joined to both B vertices and its two A
leaves joined bijectively to them. Its four columns span
$A\otimes B_V$. The remaining A vertex joined to the remaining B pair spans
$a_w\otimes B_W$. Equation (6) again gives total rank five.

In the natural four-column basis supplied by the $P_5$, the common-line vector
has every coordinate nonzero: the only non-obvious coefficient is a
$2\times2$ determinant comparing the two distinct noncentral A points, hence
is nonzero. Thus deleting any $P_5$ edge destroys the intersection; deleting a
$P_3$ edge leaves a B direction outside $B_V$. All deletions have rank five.
The graph components recover every parameter, and

\[
\boxed{N_{44}^{(6),\mathrm{ff}}=144\binom n4^2}. \tag{8}
\]

Indeed, after choosing the vertex sets, choose the $P_3$ A center (4), its B
pair (6), the $P_5$ A center (3), and the leaf-to-B bijection (2).

### 4.3 Profile $(5,4)$: two B-centered $P_3$'s plus an A-centered $P_3$

Choose five A and four B vertices. One A center is joined to a pair of B
vertices. Each of the remaining two B vertices is joined to a disjoint pair
partitioning the other four A vertices. The two B-centered $P_3$ components
span $A\otimes B_V$ for a B pair $V$; the A-centered component spans
$a_w\otimes B_W$ for the complementary pair. The same intersection and
all-nonzero-coordinate argument proves rank five and deletion-minimality.

The orientation of the three components identifies the A-centered component;
the actual two B centers distinguish the ordered assignment of the two A
pairs. Hence

\[
\boxed{N_{54}^{(6),\mathrm{ff}}=180\binom n5\binom n4}. \tag{9}
\]

The factor 180 is $5\cdot6\cdot6$: A center, its B pair, and an assignment of
the remaining four A vertices in two pairs to the two remaining B centers.

The independent constructor in `run_frontier.py` agreed as a **support set**,
not just a count, with every measured circuit in profiles $(3,4),(4,4),(5,4)$
on all four preregistered instances. This finite set equality corroborates the
proved constructions; it is not used to assert that no additional graph type
can occur for arbitrary larger point sets.

The two other measured size-six profiles are exact determinant channels:

- $(4,5)$: a four-edge component with component degrees
  $(2,2)$ on A and $(2,1,1)$ on B, plus two isolated edges;
- $(5,5)$: a $P_4$ plus three isolated edges.

Each is accepted exactly by the rank-five/all-five-deletions-rank-five
predicate of §2. Their support sets and counts vary between GF(7) and GF(13),
so neither is promoted as field-free.

## 5. One proved ambient-edge size-seven family

Choose three A and four B vertices, select one vertex on each side for an
isolated edge, and place all six edges of $K_{2,3}$ on the remaining vertices.
The $K_{2,3}$ columns are a tensor basis of $A\otimes B$, so the seven-set has
rank six. Deleting the isolated edge leaves that basis. If a basis edge is
deleted, express the isolated tensor in the chosen tensor basis: its A factor
has both coefficients nonzero and its B factor has all three coefficients
nonzero, so the coefficient at the deleted basis cell is nonzero. The isolated
tensor is outside the five-dimensional coordinate hyperplane left by the
deletion. Every six-deletion is therefore independent, and the support is a
circuit.

The two graph components recover the isolated row and column, giving

\[
\boxed{N_{34}^{(7),\mathrm{ff}}=12\binom n3\binom n4}. \tag{10}
\]

The independent constructor had exact support-set equality with every measured
profile-$(3,4)$ size-seven circuit in all four sweeps.

Other size-seven profiles are not promoted merely because some were stable
across two primes. The exact predicate remains that all seven $6\times6$
deletion determinants are nonzero.

## 6. Exhaustive finite census

The following tables are exhaustive over every subset of the indicated size.
The profile convention is `(number of A-row vertices, number of B-column
vertices)`.

### 6.1 $n=4$

Both GF(7) and GF(13) gave identical support sets.

| size | universe | profile counts | total circuits |
|---:|---:|---|---:|
| 3 | 560 | $(3,1):16$ | 16 |
| 4 | 1,820 | $(1,4):4$ | 4 |
| 5 | 4,368 | $(3,4):48$ | 48 |
| 6 | 8,008 | $(3,4):72$, $(4,4):144$ | 216 |
| 7 | 11,440 | $(3,4):48$, $(4,4):384$ | 432 |

### 6.2 $n=5$

| size | GF(7) profile counts | GF(13) profile counts |
|---:|---|---|
| 3 | $(3,1):50$ | $(3,1):50$ |
| 4 | $(1,4):25$ | $(1,4):25$ |
| 5 | $(3,4):600$, $(5,5):6$ | same totals |
| 6 | $(3,4):900$, $(4,4):3600$, $(4,5):720$, $(5,4):900$, $(5,5):360$ | $(3,4):900$, $(4,4):3600$, $(4,5):504$, $(5,4):900$, $(5,5):216$ |
| 7 | $(3,4):600$, $(3,5):5400$, $(4,4):9600$, $(4,5):23760$, $(5,4):11400$, $(5,5):12600$ | $(3,4):600$, $(3,5):5400$, $(4,4):9600$, $(4,5):26352$, $(5,4):11400$, $(5,5):15768$ |

The universes swept at $n=5$ were respectively 2,300; 12,650; 53,130;
177,100; and 480,700 subsets. Totals were 606 at size five for each prime,
6,480 versus 6,120 at size six, and 63,360 versus 69,120 at size seven.

### 6.3 Exact cross-prime support comparison at $n=5$

| size/profile | GF(7) | GF(13) | intersection | GF(7)-only | GF(13)-only |
|---|---:|---:|---:|---:|---:|
| 5 $(3,4)$ | 600 | 600 | 600 | 0 | 0 |
| 5 $(5,5)$ | 6 | 6 | 2 | 4 | 4 |
| 6 $(3,4)$ | 900 | 900 | 900 | 0 | 0 |
| 6 $(4,4)$ | 3,600 | 3,600 | 3,600 | 0 | 0 |
| 6 $(4,5)$ | 720 | 504 | 216 | 504 | 288 |
| 6 $(5,4)$ | 900 | 900 | 900 | 0 | 0 |
| 6 $(5,5)$ | 360 | 216 | 72 | 288 | 144 |
| 7 $(3,4)$ | 600 | 600 | 600 | 0 | 0 |
| 7 $(3,5)$ | 5,400 | 5,400 | 5,400 | 0 | 0 |
| 7 $(4,4)$ | 9,600 | 9,600 | 9,600 | 0 | 0 |
| 7 $(4,5)$ | 23,760 | 26,352 | 20,664 | 3,096 | 5,688 |
| 7 $(5,4)$ | 11,400 | 11,400 | 11,400 | 0 | 0 |
| 7 $(5,5)$ | 12,600 | 15,768 | 10,160 | 2,440 | 5,608 |

Equality of support sets for $(3,5),(4,4),(5,4)$ at size seven is evidence,
not a proof of a universal field-free count.

## 7. Label-invariant signatures

The machine record stores every component/degree signature. A connected
component is encoded as

\[
[e;a,b;\rho;\mathbf d_A;\mathbf d_B],
\]

where $e$ is its edge count, $a,b$ are its occupied vertex counts,
$\rho=e-a-b+1$ is cycle rank, and the degree sequences are sorted. This is a
label-invariant signature, not asserted to be a complete graph-isomorphism
invariant.

At size six, the five measured profile signatures are precisely:

| profile | component signatures |
|---|---|
| $(3,4)$ | $[4;2,2;1;(2,2);(2,2)]+[2;1,2;0;(2);(1,1)]$ |
| $(4,4)$ | $[4;3,2;0;(2,1,1);(2,2)]+[2;1,2;0;(2);(1,1)]$ |
| $(4,5)$ | $[4;2,3;0;(2,2);(2,1,1)]+2K_2$ |
| $(5,4)$ | $2[2;2,1;0;(1,1);(2)]+[2;1,2;0;(2);(1,1)]$ |
| $(5,5)$ | $[3;2,2;0;(2,1);(2,1)]+3K_2$ |

At size seven there are 17 signatures in total, grouped by profile as
$1,2,2,4,3,5$ for $(3,4),(3,5),(4,4),(4,5),(5,4),(5,5)$ respectively. The
complete exact tuples and per-prime multiplicities are retained in
`controls_results.json`; notably, multiple signatures inside the
field-dependent $(4,5)$ and $(5,5)$ profiles change multiplicity. No signature
aggregation is used as a substitute for the deletion-minor predicate.

## 8. Controls, runtime, and exact residual gaps

All 121 aggregate assertions passed. Positive controls included a known
unequal crossing and an actual seven-circuit with all six-deletions
independent. Planted rejects included a larger set containing a proper fiber,
a seven-set that was dependent but not minimal, duplicated factor columns
forcing spark two, a label-valid support with a corrupted non-tensor column,
and a block-concatenation impostor of dimension five in place of the true
six-dimensional Kronecker column. Thus both ACCEPT and REJECT paths exercised
actual rank/minor logic.

The run used exact Python integers modulo $p$, true row-major Kronecker
coordinates, no imports from a frozen directory, `sys.dont_write_bytecode=True`,
`PYTHONDONTWRITEBYTECODE=1`, `nice -n 15`, and all declared thread variables
set to one. CPU use was 40.594 seconds and wall use 42.775 seconds, below the
hard 600/900-second caps. No defect occurred and there is no broken output.

The first unresolved analytic channel is already visible at size five:

1. **Size 5, profile $(5,5)$:** convert predicate (5) into a structural
   parameterization/injective exact determinant sum, or prove no simpler
   field-free form exists.
2. **Size 6, profiles $(4,5),(5,5)$:** evaluate the rank-five/deletion-minor
   determinant sums; their support sets are demonstrably field-dependent.
3. **Size 7, profiles $(3,5),(4,4),(5,4)$:** two-prime support stability is not
   an analytic exhaustion or injective parameterization and remains OPEN.
4. **Size 7, profiles $(4,5),(5,5)$:** the deletion-determinant channels are
   field-dependent and remain OPEN.
5. For arbitrary $n$, prove or refute that the measured component signatures
   exhaust each profile; this run certifies that exhaustion only for the
   preregistered $n=4,5$ instances.

Because these obligations are part of the gate's promotion rule, the single
whole-gate verdict is **FROZEN-INCONCLUSIVE**, while the complete two-row
capstone and the explicitly proved unequal-row families (3)--(10) remain valid
certified subresults.
