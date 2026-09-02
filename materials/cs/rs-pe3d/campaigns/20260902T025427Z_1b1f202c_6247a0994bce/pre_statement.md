# PREREG-RSPE3D-12-3ROW-DETERMINANTAL — gate **H-3ROW-DETERMINANTAL**

Written 2026-09-02 UTC by agent `RsPe3dH2`, after run
`20260902T023120Z_8ac826db_96f766bfe6cc` was frozen, closed
FROZEN-INCONCLUSIVE, independently verified by Main, and path-scoped committed.
This file is committed before any new parameter-dependent computation. Binding
order: preregistration -> path-scoped commit -> `campaign.py init` ->
byte-identical in-run copy with source commit and SHA-256 -> controls/censuses
-> proof record -> `freeze` -> exactly one `close --verdict`. Write only under
`cs/rs-pe3d/`; never edit a frozen run or root ledger.

The candidate route below was analytically derived from the bi-Vandermonde
matrix and the frozen predecessor census before this preregistration, but it is
not authority. Every set equality and formula is an obligation. A failed
reduced predicate is preserved as a defect/counterexample and is not promoted.

## 1. Setting and exact matrices

Let

\[
a_x=(1,x)^{\mathsf T},\qquad b_y=(1,y,y^2)^{\mathsf T},\qquad
h(x,y)=(1,y,y^2,x,xy,xy^2)^{\mathsf T},
\]

for distinct evaluation points within finite sets $X,Y\subset F$. Thus
$(d_A,d_B)=(3,4)$ and the product ambient dimension is six. A indices are row
vertices and B indices are column vertices, so every nonfiber circuit of size
5--7 has

\[
|\pi_A S|\ge3,
\quad |\pi_BS|\ge4,
\quad \deg_A\le3,
\quad \deg_B\le2.
\]

For an ordered six-set $T=\{(x_i,y_i):1\le i\le6\}$ define

\[
D_6(T)=\det[h(x_i,y_i)]_{i=1}^6.
\]

Laplace expansion along the first three rows gives the explicit polynomial

\[
D_6(T)=\sum_{I\subset[6],\ |I|=3}(-1)^{\sum_{i\in I}i}
\left(\prod_{j\notin I}x_j\right)V(y_I)V(y_{I^c}),
\tag{L}
\]

where $V$ is the ordered $3\times3$ Vandermonde determinant. For four points
write

\[
D_4(R)=\det[1,y,x,xy]_{(x,y)\in R}.
\]

All machine decisions use exact integer arithmetic modulo $p$ and the actual
six-coordinate columns. Formula (L) and reduced predicates are evaluated
independently of direct rank where set equality is claimed.

## 2. P1 — complete size-five classification

Prove that a size-five circuit is exactly one of:

1. a profile-$(3,4)$ Theorem-X crossing, already field-free, count
   $12\binom{|X|}{3}\binom{|Y|}{4}$;
2. a profile-$(5,5)$ matching lying on the graph of a unique nondegenerate
   Möbius map $M:y\mapsto x$.

For five distinct B coordinates, use barycentric weights

\[
w_j=\left(\prod_{k\ne j}(y_j-y_k)\right)^{-1}
\]

and prove

\[
\ker[1;y;y^2]=\operatorname{span}\{(w_j),(w_jy_j)\}.
\]

A product relation must have

\[
\lambda_j=w_j(\alpha+\beta y_j),
\qquad
x_j(\alpha+\beta y_j)=\gamma+\delta y_j,
\]

which is exactly the Möbius graph condition; circuit minimality excludes a
pole/degenerate constant. With

\[
K_M=\{y\in Y:M(y)\in X\},\qquad k_M=|K_M|,
\]

prove the injective exact sum

\[
\boxed{N_{55}^{(5)}(X,Y)=\sum_{M\in\mathrm{PGL}(2,F)}\binom{k_M}{5}}.
\tag{P1}
\]

For four distinct B coordinates with one repeated twice, use the two-dimensional
B-kernel to prove that a full-support five-relation forces the three singleton
B edges to share one A vertex, hence is precisely a crossing. This must exclude
all other profiles analytically, not by census.

## 3. P2 — six-set determinant factorization and residual templates

Prove the square-evaluation lemma. If a six-set has one B value $y_0$ repeated
at two distinct A values $x_a,x_b$, and four singleton-B edges $R$, then, up to
the sign fixed by column order,

\[
D_6(S)=\pm(x_b-x_a)\prod_{(x,y)\in R}(y-y_0)D_4(R).
\tag{F}
\]

Equivalently, a nonzero polynomial $f(y)+xg(y)$ of bidegree at most $(1,2)$
vanishing on the repeated pair must have both $f,g$ divisible by $y-y_0$;
the quotient is the bilinear/Möbius condition on $R$. Verify (F) numerically
as an exact identity, not only its zero set.

Use P1 to prove both residual circuit predicates.

### P2a. Profile $(4,5)$

The graph must be $P_5+2K_2$. Its four singleton-B edges form an all-distinct
matching $R$. Prove

\[
S\text{ is a circuit}\iff D_4(R)=0.
\tag{P2a}
\]

All five-deletions are structurally independent. Parameterize $R$ by its unique
$M$, choose a new B value outside its four-domain set, and choose two of the
four matching edges to receive the new common B neighbor. Prove injectivity and

\[
\boxed{N_{45}^{(6)}(X,Y)=
6(|Y|-4)\sum_M\binom{k_M}{4}}.
\tag{C45}
\]

### P2b. Profile $(5,5)$

The graph must be $P_4+3K_2$. Let $R$ be the four matching edges left after
removing both edges incident to the unique degree-two B vertex; let $Q$ be the
five-edge matching obtained by deleting the central edge joining the unique
degree-two A and B vertices. Prove

\[
S\text{ is a circuit}\iff D_4(R)=0\ \text{and}\ Q
\text{ is not a size-five Möbius circuit}.
\tag{P2b}
\]

The second clause is the only nonautomatic proper-subset condition. Recover all
parameters from the graph and prove

\[
\boxed{N_{55}^{(6)}(X,Y)=4\sum_M\binom{k_M}{4}
\bigl((|X|-4)(|Y|-4)-(k_M-4)\bigr)}.
\tag{C55}
\]

Also prove the unified field-free $q=4$ size-six count. If the B-degree pattern
is $(2,2,1,1)$, dependence is equivalent to the two singleton-B edges sharing
an A row; minimality is equivalent to that row avoiding both doubled-B pairs.
For profile $(p,4)$, $p=3,4,5$, put

\[
h_q=\sum_{t=0}^q(-1)^t\binom qt\binom{q-t}{2}^2.
\]

Then

\[
N_{p4}^{(6)}=6p\binom{|X|}{p}\binom{|Y|}{4}h_{p-1},
\tag{Q4-6}
\]

recovering coefficients 18, 144, 180.

## 4. P3 — explicit size-seven deletion-minor predicates

A seven-set is automatically dependent. Prove it is a circuit iff every
six-deletion has nonzero determinant, and reduce each determinant according to
its B-degree profile. The following is the preregistered complete predicate.

1. **$q=4$, B degrees $(2,2,2,1)$.** Deleting the unique singleton-B edge
   leaves three B basis columns each paired with two distinct A columns and is
   automatically independent. Deleting an edge at a doubled B vertex reduces
   to the $q=4$ six-set lemma. All minors are nonzero iff the A endpoint of the
   unique singleton-B edge occurs in none of the three doubled-B pairs.
   For profile $(p,4)$, $3\le p\le7$, put
   \[
   g_q=\sum_{t=0}^q(-1)^t\binom qt\binom{q-t}{2}^3.
   \]
   Prove the injective field-free count
   \[
   \boxed{N_{p4}^{(7)}=4p\binom{|X|}{p}\binom{|Y|}{4}g_{p-1}}.
   \tag{Q4-7}
   \]
2. **$q=5$, B degrees $(2,2,1,1,1)$.** Let the three singleton-B edges have
   A endpoints $r_1,r_2,r_3$. Require them pairwise distinct. For each endpoint
   of each doubled-B pair, form the four-matching consisting of that endpoint
   at its B value together with the three singleton-B edges; require all four
   corresponding $D_4$ values nonzero. These three combinatorial and four
   minor conditions are necessary and sufficient.
3. **$q=6$, B degrees $(2,1,1,1,1,1)$.** For each singleton-B deletion,
   require the $D_4$ of the four remaining singleton-B edges to be nonzero.
   For each of the two doubled-B-edge deletions, require the all-distinct-B
   $D_6$ from (L) to be nonzero.
4. **$q=7$, all B degrees one.** Require all seven all-distinct-B deletion
   determinants (L) to be nonzero.

For all-distinct B six-sets, record the equivalent condition: $D_6=0$ iff the
six paired points lie on a possibly degenerate rational graph
$x=-f(y)/g(y)$ with $\deg f,\deg g\le2$. Do not call this a unique
parameterization without proving uniqueness/nondegeneracy.

For every measured size-seven profile and label-invariant component/degree
signature, record direct and predicted support counts. Overall and per-signature
set equality is required. Prime-stable counts alone never prove a formula.

## 5. General circuit principle, exact scope

Prove only the linear-matroid statement actually justified. For an
$m$-column support $S$ in rank at most six,

\[
S\text{ is a circuit}\iff \operatorname{rank}H_S=m-1
\text{ and every one-column deletion has rank }m-1.
\]

For $m=5,6,7$ this becomes respectively: all $5\times5$ maximal minors vanish
plus nonzero deletion minors; $D_6=0$ plus nonzero five-deletion minors; and
seven nonzero $D_6$ deletion determinants. This is a minimal rank-drop
criterion for the actual Kronecker evaluation matrix. Do not promote a purely
graphical version when field-dependent minors remain.

## 6. Exhaustive matrix and controls

After initialization, independently census $X=Y=\{1,\ldots,n\}$ for
$n\in\{4,5\}$ and $p\in\{7,11,13\}$ at sizes 5, 6, 7. Use true row-major
Kronecker columns. For every configuration:

- identify direct circuits by exact rank and every deletion;
- independently construct/predict P1, P2, and P3 supports without consulting
  the direct circuit set;
- assert exact support-set equality overall and by profile/signature;
- assert (L), (F), PGL sums, and all formulas where applicable;
- compare prime support sets exactly and retain intersections/one-sided
  differences.

Controls in both directions:

1. known crossing and Möbius-matching ACCEPTs;
2. an inserted proper crossing or Möbius five-circuit inside a larger set must
   REJECT;
3. the deliberately wrong profile-$(5,5)$ size-six predicate `D4(R)=0` without
   the central-deletion condition must overaccept, and exact set equality must
   fail; plant the identity-Möbius example explicitly;
4. dependence of an arbitrary seven-set is automatic but minimality is tested
   separately; include one actual seven-circuit ACCEPT and one nonminimal
   REJECT;
5. true six-coordinate Kronecker versus five-coordinate block-concatenation
   guard;
6. a label-valid support with one corrupted actual product column must fail the
   appropriate actual minor predicate.

Any mismatch writes and preserves `broken_results.json`, discloses the failing
family, and reruns every affected family after a justified fix. No broken
artifact is overwritten.

## 7. Runtime, promotion, and lifecycle

Use Python stdlib only, exact integers modulo $p$,
`sys.dont_write_bytecode=True`, `PYTHONDONTWRITEBYTECODE=1`, `nice -n 15`, and
`OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=VECLIB_MAXIMUM_THREADS=1`.
Hard budget: CPU 900 seconds, wall 1200 seconds. Never import from a frozen run.

Close **FROZEN-CERTIFIED** only if P1, P2a, P2b, factorization (F), the complete
P3 deletion-minor predicates, and all six direct-versus-predicted set equalities
are proved with no gap. A closed determinant sum need not be simplified to a
field-free polynomial, but every claimed injective sum must be proved. If any
channel lacks a proved explicit predicate or any set mismatch survives the
defect pass, close **FROZEN-INCONCLUSIVE** and name the first exact gap and the
strongest exhaustive range. Freeze before exactly one close. After closing,
update only `SCOPE_NOTE_H_MIX_CROSSING.md` and `state.json`, then path-scoped
commit; root ledgers remain Main's responsibility.
