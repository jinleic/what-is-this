# PREREG-RSPE3D-10-DP2-SIZE5 — gate **H-DP2-SIZE5**

Written 2026-09-01 (local date), agent `RsPe3dH2`, after
H-DP1-THREEFACTOR was frozen, closed FROZEN-CERTIFIED, independently verified
by Main, recorded as the root-ledger capstone, and committed. This file is
path-scoped committed before any parameter-dependent compute for this run.
Binding order: preregistration → `git commit --only` → `campaign.py init
--gate H-DP2-SIZE5` → byte-identical prereg copy in the minted run directory
with source commit and SHA-256 → exact controls and censuses → proof record →
`freeze` → `close --verdict` with exactly one verdict. No frozen campaign
directory or root ledger is modified.

The candidate formulas, decomposition description, and numeric anchors below
were supplied by Main as derived and externally verified, not as authority.
This run must prove the support classification, independently reconstruct the
families, reproduce every anchor, and reject any oversimplified candidate law.

## 1. Fixed setting and target

Let $A,B$ be $2\times n$ GRS/Vandermonde parity-check factors over $F_p$, with
columns

$$a_x=(1,x)^{\mathsf T},\quad x\in X,
\qquad b_y=(1,y)^{\mathsf T},\quad y\in Y,$$

where $X,Y$ contain $n$ distinct field points. Both sparks are 3. Product
columns are

$$h(x,y)=a_x\otimes b_y=(1,y,x,xy)^{\mathsf T}.$$

A five-cell set $S\subseteq X\times Y$ is a circuit iff its five product
columns have rank four and every four-subset has rank four. Its profile is
$(|\pi_X(S)|,|\pi_Y(S)|)$.

The target is a complete size-five support theorem:

1. the six profiles with one projection of size three are exactly the tight
   $(1,2)$ or $(2,1)$ cancellation families and have field-free closed counts;
2. the remaining profiles $(4,4),(4,5),(5,4),(5,5)$ have complete graph
   templates and exact determinant predicates (plus a closed reduction of the
   $(4,4)$ count to the frozen size-four all-distinct count);
3. every possible profile/template is either accepted by those predicates or
   rejected because it contains a proper size-three fiber or size-four
   crossing/all-distinct circuit.

## 2. Correct support-cardinality law and decomposition obligations

Write a product relation as $\Gamma=C_1+C_2$, with every nonzero column of
$C_1$ in $\ker A$ and every nonzero row of $C_2$ in $\ker B$. Put

$$P=\operatorname{supp}C_1,
\quad Q=\operatorname{supp}C_2,
\quad r=|P\cap Q|,
\quad c=|\{z\in P\cap Q:C_1(z)+C_2(z)=0\}|.$$

Let $T_1=|P|$ and $T_2=|Q|$. The exact general law to prove is

$$\boxed{|\operatorname{supp}\Gamma|=T_1+T_2-r-c.}\tag{L}$$

Indeed every overlap has been counted twice in $T_1+T_2$ and must be counted
once if it survives, zero times if it cancels. The supplied law

$$|S|=k_1d_A+k_2d_B-2c\tag{L-tight}$$

is a **corollary**, not the unrestricted law: it requires exactly $k_1,k_2$
minimum-support relations ($T_1=k_1d_A$, $T_2=k_2d_B$) and every overlap to
cancel ($r=c$).

For $d_A=d_B=3$, call such a tight/full-cancel decomposition irreducible if
every active row and column contributes at least one surviving cell. Then each
active minimum relation has at most two canceled overlaps, so
$c\le2k_1$ and $c\le2k_2$. Combine this with

$$5=3(k_1+k_2)-2c$$

to prove that the only irreducible solutions are

$$(k_1,k_2,c)=(1,2,2)\quad\text{or}\quad(2,1,2).$$

Larger formal integer solutions are padded/reducible (some active minimum
relation has no survivor) and are not circuit shapes. Remaining five-circuits
must use a noncanceled overlap and/or a relation of support greater than three.
The theorem must state this correction explicitly; it must not promote
(L-tight) as a general identity.

## 3. Complete graph/minor classification

A five-circuit contains no fiber triple, so every row/column degree is at most
two. With five edges, every projection therefore has size at least three. Use
the frozen complete size-four theorem: once fiber triples are excluded, a
four-subset is dependent iff it is a profile-$(3,3)$ crossing or an
all-distinct profile-$(4,4)$ set with zero bilinear determinant. Prove:

> A five-cell support is a circuit iff (i) all row/column degrees are at most
> two, (ii) no four-edge subset is a crossing, and (iii) every all-distinct
> four-edge subset has nonzero determinant.

Enumerate every degree-at-most-two bipartite graph template (paths use edge
length):

- $(3,3)$: $C_4+K_2$ (accepted) or $P_6$ (contains a crossing);
- $(3,4)$: $P_5+K_2$ (accepted) or $P_4+P_3$ (contains a crossing);
- $(3,5)$: two row-centered $P_3$ components plus $K_2$ (accepted);
- mirrors $(4,3),(5,3)$;
- $(4,4)$: $P_4+2K_2$ (one determinant) or two oppositely centered
  $P_3$ components plus $K_2$ (contains a crossing);
- $(4,5)$: one row-centered $P_3+3K_2$ (two determinants), and its $(5,4)$
  mirror;
- $(5,5)$: $5K_2$ (five determinants).

This list must be proved exhaustive from component/degree counts, not inferred
from census.

## 4. The tight $(1,2)$ family and its mirror

Choose:

- a row 3-circuit $R=\{u_1,u_2,w\}$;
- one active $C_1$ column $j_0$ supported on $R$;
- active $C_2$ rows $u_1,u_2$;
- row 3-circuits $Z_1,Z_2$, each containing $j_0$;
- scales that cancel exactly $(u_1,j_0),(u_2,j_0)$.

The support is

$$S=\{(w,j_0)\}
\cup\bigl(\{u_1\}\times(Z_1\setminus\{j_0\})\bigr)
\cup\bigl(\{u_2\}\times(Z_2\setminus\{j_0\})\bigr).\tag{F}$$

Prove directly from the templates that (F) is a genuine circuit: no row or
column has degree three, it contains no four-edge crossing, and it has no
all-distinct four-subset. In particular it has no size-three fiber and no
size-four crossing. Build exact relation coefficients from the unique factor
triple relations and verify both cancellations.

Prove parameter injectivity by recovering the unique isolated edge
$(w,j_0)$, the two active rows, and their incident column sets. Derive:

1. $Z_1=Z_2$ (profile $(3,3)$):
   $$N_{33}=9\binom n3^2.$$
2. $|Z_1\cap Z_2|=2$ (profile $(3,4)$):
   $$N_{34}=18n\binom n3\binom{n-1}3
   =72\binom n3\binom n4.$$
3. $Z_1\cap Z_2=\{j_0\}$ (profile $(3,5)$):
   $$N_{35}=18n\binom n3\binom{n-1}4.$$

The $(2,1)$ mirror gives $N_{33}$, $N_{43}=N_{34}$, and $N_{53}=N_{35}$.
The same profile-$(3,3)$ support is generated by both orientations but counted
once as a set; prove the two descriptions coincide. Explain the numerical
coincidence $N_{33}=9\binom n3^2$ with the size-four crossing count: both
families choose the same $(R,Z,w,j_0)$ parameters, but size four uses the two
cross arms while size five uses the complementary $2\times2$ rectangle plus
the opposite isolated edge.

These counts are field-independent because only distinctness/spark three and
combinatorial choices enter.

## 5. Remaining profiles and canonical decomposition shapes

### P44 — profile $(4,4)$

A surviving support is $P_4+2K_2$. Its unique edge joining the degree-two row
to the degree-two column is the central edge. Removing it leaves one
all-distinct four-edge matching $Q$. Prove

$$S\text{ is a circuit}\iff \det H_Q\ne0.$$

Conversely every independent all-distinct four-matching $Q$ and every ordered
pair of distinct edges of $Q$ yield one unique $S$ by adding the cross edge
from the row of the second edge to the column of the first. Therefore, with

$$M_4=24\binom{|X|}4\binom{|Y|}4$$

and $N_{\rm all}^{(4)}(X,Y)$ the frozen size-four all-distinct circuit count,
prove the closed reduction

$$\boxed{N_{44}^{(5)}=12\bigl(M_4-N_{\rm all}^{(4)}(X,Y)\bigr).}$$

Give a canonical $(2,2)$ decomposition: two minimum-support $A$ relations and
two minimum-support $B$ relations meet in four cells; three overlap cells
cancel and the central overlap survives. Thus
$(T_1,T_2,r,c)=(6,6,4,3)$ and (L), not (L-tight), gives five.

### P45 — profiles $(4,5)$ and $(5,4)$

A $(4,5)$ candidate is a row-centered $P_3+3K_2$. It has exactly two
all-distinct four-subsets, obtained by retaining one or the other edge at the
degree-two row. Prove it is a circuit iff both determinants are nonzero.
Equivalently sum this two-minor predicate over: four chosen rows, five chosen
columns, the repeated row, its unordered two columns, and a bijection from the
other three rows to the remaining three columns. State the mirror predicate
for $(5,4)$.

Give a canonical $(2,2)$ decomposition with all four overlaps canceled and
one side's total relation support seven rather than six:
$(T_1,T_2,r,c)=(6,7,4,4)$ for $(4,5)$ and $(7,6,4,4)$ for $(5,4)$. Construct
it by choosing the second active row/column, using three minimum factor
relations, and solving the remaining two coefficients against an independent
factor pair. Show circuit minimality forces every intended survivor coefficient
nonzero.

### P55 — profile $(5,5)$

A candidate is a five-edge matching. Prove it is a circuit iff all five
four-edge matching minors are nonzero, equivalently the five product columns
have rank four and their unique relation has full support. The factor-rank
inequality $\rho_A+\rho_B\le5$ is necessary but **not sufficient** here: for
2-row factors it accepts every matching ($2+2\le5$), including matchings with
a zero four-minor. This insufficiency must be planted and rejected.

Give a canonical decomposition of every accepted matching: choose two matched
$A$ columns as a basis and apply their dual coordinate functionals to the full
product relation. This creates two $B$-kernel rows of support four; the other
three matched columns become three minimum-support $A$-kernel columns. Six
overlap cells cancel, giving
$(k_1,k_2,T_1,T_2,r,c)=(3,2,9,8,6,6)$ and
$9+8-6-6=5$ (or the mirror $(2,3)$ construction).

For P45/P55, an exact finite determinant sum is a complete count; no
field-free polynomial closed form is claimed.

**Promotion rule.** FROZEN-CERTIFIED only if (L), the tight-shape enumeration,
the graph/minor theorem, all six closed formulas/injectivity proofs, P44's
closed reduction, P45/P55 predicates, and the stated canonical decompositions
are proved without a gap; all exhaustive support sets agree exactly at every
registered prime; and every plant fires. Otherwise close FROZEN-INCONCLUSIVE
with the first exact gap and strongest exhaustive range.

## 6. Exact control matrix

All mathematical decisions use exact Python integers modulo $p$, actual
$(1,y,x,xy)$ product columns, exact Gauss–Jordan rank/determinants, and stdlib
only. A direct circuit route checks rank four and every four-subset; an
independent structural route checks graph templates, forbidden crossings, and
only the prescribed all-distinct minors. Full support-set equality is required.

### T1 — $n=4$, three primes

For $X=Y=\{1,2,3,4\}$ and $p\in\{7,11,13\}$ exhaust all
$\binom{16}{5}=4368$ five-subsets. Pin supplied anchors:

- GF(13): 864 circuits with profiles
  $(3,3):144$, $(3,4):288$, $(4,3):288$, $(4,4):144$;
- GF(7): 912 with the same first three counts and $(4,4):192$.

At GF(11), record the exact determinantal count without a pre-supplied pin.
At every prime assert direct/structural set equality, exact equality of the
six field-free constructed families with all measured profiles having a
projection of size three, and
$N_{44}=12(M_4-N_{\rm all}^{(4)})$ as sets and counts.

### T2 — $n=5$, three primes

For $X=Y=\{1,2,3,4,5\}$ and $p\in\{7,11,13\}$ exhaust all
$\binom{25}{5}=53130$ five-subsets. Pin GF(13): total 17880 with

$$(3,3):900,
(3,4):3600,
(3,5):900,
(4,3):3600,
(4,4):6192,
(4,5):864,
(5,3):900,
(5,4):864,
(5,5):60.$$

At GF(7)/GF(11), record exact determinantal channels without pre-supplied
counts. At every prime assert:

- direct/structural full support-set equality;
- field-free profile counts 900/3600/900 and mirrors;
- exact set equality with constructed $(1,2)/(2,1)$ families;
- P44 construction/set equality and count reduction;
- P45/P54 two-minor constructor set equality;
- P55 five-minor matching constructor set equality;
- canonical decomposition reconstruction on every surviving P44/P45/P54/P55
  support.

### T3 — controls in both directions

1. **Known-true ACCEPT.** Materialize one support from each of the tight
   $(1,2)$ profiles $(3,3),(3,4),(3,5)$ at $n=5$, verify exact decomposition
   coefficients, rank four, and all five deletion minors nonzero.
2. **Inserted-crossing REJECT.** Take a known size-four crossing and add a
   fifth distinct cell. Assert the five-set contains that dependent proper
   subset and is rejected by both rank and structural routes.
3. **Factor-rank-only REJECT.** At $n=5,p=13$, assert all 120 all-distinct
   matchings satisfy $\rho_A+\rho_B=4\le5$ but only the exact determinant
   subset (pinned 60) are circuits.
4. **Non-tensor label REJECT.** On a label-valid tight support, replace one
   true product column by a duplicate of another actual column. The label
   predicate still accepts; actual rank/proper-subset checks reject.
5. **Duplicated-factor-column REJECT.** Clone one factor column; assert spark
   drops from 3 to 2 and a dependent pair appears.
6. **Concat-vs-Kronecker guard.** For 2-row by 2-row factors assert true
   dimension four and exact multiplicative coordinates; planted block
   concatenation has the wrong coordinates even though its dimension is also
   four.
7. **Wrong general-law REJECT.** Exhibit a certified P44 canonical
   decomposition with $(T_1,T_2,r,c)=(6,6,4,3)$: (L) gives 5 while replacing
   $r$ by $c$ and using the unrestricted `-2c` formula gives 6. Assert the
   latter is rejected.

## 7. Arithmetic, budget, defects, and kill criteria

- Exact field arithmetic only; inverses via `pow(x,-1,p)`; no floating-point
  mathematical decisions, NumPy, or external algebra. Timing does not enter a
  theorem decision.
- `sys.dont_write_bytecode=True`; no frozen import; execute with
  `PYTHONDONTWRITEBYTECODE=1`.
- One process at `nice -n 15`, with
  `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=VECLIB_MAXIMUM_THREADS=1`.
- Hard budget: CPU 600 seconds, wall 900 seconds. Abort rather than raise it.
- Every positive claim has a known-true ACCEPT and a planted REJECT in-run.
- Any assertion mismatch is fail-loud and writes `broken_results.json` before
  nonzero exit. Preserve every broken artifact. Disclose defects, repair the
  cause without weakening pins, and rerun every affected family.
- No parameter-dependent compute may occur before this file's path-scoped
  commit. If amended after init, provenance records both hashes.
- The graph/minor proof is mandatory: census cannot promote an unproved
  decomposition taxonomy. Conversely, no field-free formula may be claimed
  for a determinantal channel merely because three primes happen to agree.
