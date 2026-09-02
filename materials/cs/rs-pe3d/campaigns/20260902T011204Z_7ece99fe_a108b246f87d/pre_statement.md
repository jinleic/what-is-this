# PREREG-RSPE3D-7-ALLDISTINCT-CROSSRATIO — gate **H-ALLDISTINCT-CROSSRATIO**

Written 2026-09-01 (local date), agent `RsPe3dH2`, after the sequential
predecessor H-MIX-H2-CLOSURE was frozen, closed FROZEN-CERTIFIED, and committed.
This file is path-scoped committed before any parameter-dependent compute for
this run. Binding order: preregistration → `git commit --only` of this file →
`campaign.py init --gate H-ALLDISTINCT-CROSSRATIO` → byte-identical copy into
the minted run directory with source commit and SHA-256 → exact controls and
censuses → proof record → `freeze` → `close --verdict` with exactly one verdict.
No frozen campaign directory is modified.

The characterization and numeric anchors below were supplied by Main as a
derived-and-externally-verified candidate, not as authority. This run must
verify every algebraic implication independently. A failed implication leaves
the all-distinct characterization open regardless of numerical agreement.

## 1. Fixed setting and statement to prove

Let $A$ be the $2\times n_A$ matrix with columns $a_u=(1,x_u)^{\mathsf T}$
and $B$ the $2\times n_B$ matrix with columns $b_v=(1,y_v)^{\mathsf T}$ over a
field, with the $x_u$ pairwise distinct and the $y_v$ pairwise distinct.
Thus $d_A=d_B=3$, $d=3$, and a product column is the row-major pure tensor

$$h(u,v)=a_u\otimes b_v=(1,y_v,x_u,x_u y_v)^{\mathsf T}\in F^4.$$

An **all-distinct support** of size four is
$S=\{(u_k,v_k):k=1,2,3,4\}$ with all $u_k$ distinct and all $v_k$ distinct.
Fix an ordering of the four pairings and define

$$\operatorname{CR}(a,b,c,d)=\frac{(a-c)(b-d)}{(a-d)(b-c)},$$

whose denominators are nonzero on distinct quadruples.

**Theorem AD-CR (fixed target).** The following are equivalent for an
all-distinct size-four support:

1. $S$ is a circuit of $A\otimes B$.
2. $\det[(1,y_k,x_k,x_k y_k)]_{k=1}^4=0$ (the four product columns are the
   columns of the displayed $4\times4$ matrix).
3. A nonzero bilinear form
   $c_0+c_1y+c_2x+c_3xy$ vanishes at all four paired points.
4. The paired points are the graph of one nondegenerate Möbius transformation
   $y=-(c_2x+c_0)/(c_3x+c_1)$, with no pole at a selected $x_k$.
5. The simultaneous pairing preserves cross-ratio:
   $$\operatorname{CR}(x_{u_1},x_{u_2},x_{u_3},x_{u_4})=
   \operatorname{CR}(y_{v_1},y_{v_2},y_{v_3},y_{v_4}).$$

Equality is invariant under a simultaneous reordering of the four pairs, so
item 5 is a support property. With $u_1<\cdots<u_4$ fixed for counting,

$$N_{\rm alldistinct}=
\#\{(U,\sigma): |U|=4,\ \sigma:U\hookrightarrow[n_B],\
\operatorname{CR}(x_U)=\operatorname{CR}(y_{\sigma(U)})\},$$

where an injection is the pairing and hence counts each support once. The
field dependence of the channel is exactly the field dependence of these
cross-ratio coincidences for the reduced point configurations.

## 2. Proof obligations (all must appear in `theorem_alldistinct_crossratio.md`)

- **P1 — product columns and circuit minimality.** Derive
  $(1,y,x,xy)$ from the tensor convention. Prove any three product columns
  with pairwise-distinct $x$ and pairwise-distinct $y$ are independent (for
  example, annihilate one $A$ factor and use independence of the other two
  $B$ factors). Therefore determinant zero for four all-distinct columns is
  dependence of rank exactly three, hence a genuine circuit, not merely a
  dependent four-set.
- **P2 — determinant iff bilinear vanishing.** Prove a square $4\times4$
  matrix has zero determinant iff it has a nonzero left-null vector
  $(c_0,c_1,c_2,c_3)$, which is exactly simultaneous vanishing of
  $c_0+c_1y+c_2x+c_3xy$.
- **P3 — bilinear curve iff Möbius graph, including degeneracies.** Write
  $D(x)=c_3x+c_1$, $N(x)=c_2x+c_0$ and
  $D(x)y+N(x)=0$. Prove $D\not\equiv0$; prove no selected $x_k$ is a pole;
  and prove $\Delta=c_0c_3-c_1c_2\ne0$. Explicitly handle
  $c_1c_2=c_0c_3$ (the degenerate/factored case) and show it contradicts four
  distinct $x$'s and $y$'s. If $c_3=0$, show the curve is a non-axis-parallel
  affine line. If $c_3\ne0$, show
  $$(x+c_1/c_3)(y+c_2/c_3)=
  (c_1c_2-c_0c_3)/c_3^2=:r,$$
  with $r\ne0$, and show a selected point cannot sit at the pole
  $x=-c_1/c_3$. Conversely translate every nondegenerate Möbius graph back to
  a bilinear form.
- **P4 — Möbius graph iff cross-ratio equality.** Prove invariance using
  $M(s)-M(t)=\det(M)(s-t)/((cs+d)(ct+d))$. Conversely construct the unique
  Möbius transformation sending three distinct $x$'s to three distinct
  $y$'s; use cross-ratio equality and injectivity in the fourth argument to
  force the fourth pair. Address simultaneous reordering and poles.
- **P5 — count and exact scope.** Prove the injection formula counts each
  all-distinct support once and that field variation is precisely variation
  in cross-ratio coincidences. For general $r_A\times n_A$, $r_B\times n_B$
  GRS factors, prove only the safe statement: product columns, up to nonzero
  GRS scalars, are evaluations of the bi-Vandermonde monomials
  $\{x^iy^j:0\le i<r_A,0\le j<r_B\}$; dependence is a column-rank condition on
  that matrix. When $r_Ar_B>d+1$ it is not a single square determinant, and
  circuitness additionally requires every proper selected submatrix to have
  full column rank. No general higher-rank cross-ratio formula is claimed.

**Promotion rule.** FROZEN-CERTIFIED only if P1–P5 are each discharged with no
gap (especially the degenerate determinant-zero and pole cases), all exact
set-equality controls pass, and every planted REJECT fires. Otherwise close
FROZEN-INCONCLUSIVE, naming the first exact proof or control gap.

## 3. Exact control matrix

All controls run after `init`, with Python integers modulo $p$, exact
Gauss–Jordan rank/determinant/nullspace, stdlib only.

### T1 — symmetric 4×4 point pair, five primes

For $x=y=(1,2,3,4)$ and each $p\in\{7,11,13,17,31\}$:

1. Enumerate all 24 bijections between the four $x$ and four $y$ points.
2. Assert support-level equality among (i) determinant-zero pairings,
   (ii) cross-ratio-equal pairings, and (iii) the full product census's
   profile-$(4,4)$ circuit population.
3. Pin the all-distinct counts to $8,4,12,4,4$, respectively, and pin the full
   size-four circuit totals to $152,148,156,148,148$, with profile-$(3,3)$
   crossing count 144 at every prime. In particular GF(31) must give
   $148=144+4$ exactly.
4. On every determinant-zero pairing, solve a nonzero bilinear left-null
   vector and assert exact vanishing, nonzero Möbius determinant, and no pole.

The identity pairing and reversal are known-true pure-tensor ACCEPTs at every
prime (lines $y=x$ and $x+y=5$). At GF(13), additionally assert
$\{(1,1),(2,3),(3,4),(4,2)\}$ is a genuine hyperbola circuit, not an affine-line
case; extract and verify its nonzero $r$.

### T2 — asymmetric 5×5 pair, two primes

For $x=(1,2,3,4,5)$, $y=(2,3,5,7,11)$ at $p=13$ and $p=17$, enumerate every
all-distinct support: $\binom54\binom54 4!=600$ pairings per prime. Assert
set equality of determinant-zero, cross-ratio-equal, and exhaustive-census
profile-$(4,4)$ populations; assert zero disagreement across all 1,200
pairings; record (do not pre-pin) the resulting counts. Check the bilinear /
Möbius witness conditions on every accepted support.

### T3 — planted REJECTs and tensor-layout guard

1. **Cross-ratio mismatch REJECT (GF(13)).** Pair $x=(1,2,3,4)$ with
   $y=(1,2,4,3)$; assert the two cross-ratios differ and the determinant is
   nonzero, so the support is rejected.
2. **Non-tensor plant REJECT (GF(13)).** In the symmetric product replace
   actual column $h(0,0)$ by $h(0,0)+h(0,1)$ while retaining the original grid
   labels. The pristine identity pairing is a CR-equal circuit; the planted
   actual columns must become independent, creating a fail-loud disagreement
   with the label-based criterion and rejecting the non-tensor instrument.
3. **Duplicated-column REJECT.** Clone one factor column; assert factor spark
   drops from 3 to 2 and a dependent pair is detected.
4. **Concat-vs-Kronecker guard.** On a 2-row by 3-row factor pair, assert the
   independently formed row-major Kronecker columns have dimension 6 and exact
   multiplicative coordinates, while planted block concatenation has dimension
   5 and different coordinates.

## 4. Exact arithmetic, budgets, defects, and kill criteria

- Exact mathematical decisions only: Python integers modulo $p$; inverse by
  `pow(x,-1,p)`; no floats, NumPy, or external algebra package. Timing values
  do not enter mathematical decisions.
- `sys.dont_write_bytecode=True` before any possible frozen-directory import;
  the run imports nothing from frozen directories. Execute with
  `PYTHONDONTWRITEBYTECODE=1`.
- One process under `nice -n 15`, with
  `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=VECLIB_MAXIMUM_THREADS=1`.
- Hard cap: 600 CPU seconds and 900 wall seconds, checked throughout each
  exhaustive family; no adaptive extension beyond the fixed 5-prime and
  2-prime matrices above.
- Any defect is disclosed in-run, broken outputs preserved, and every affected
  family rerun from its beginning. No failed assertion may be weakened or
  silently repinned.
- Kill on any P1–P5 gap; determinant/CR/census set disagreement; anchor
  mismatch; ACCEPT failure; REJECT failure; incomplete enumeration; tensor
  layout, bytecode, thread, CPU, or wall-budget violation; or unresolved
  defect. Close exactly one verdict.

## 5. Scope and nonclaims fixed in advance

Certified if promoted: the 2-row GRS all-distinct size-four characterization,
its exact cross-ratio count formula, and field-dependence mechanism, with exact
finite corroboration on the listed configurations. Not certified: a closed
formula evaluating the cross-ratio coincidence count for arbitrary point sets;
spark-2 branches; a cross-ratio criterion for $r_Ar_B>4$; general
$r_A,r_B$ circuit classification beyond the safe bi-Vandermonde rank statement;
three or more factors; repeated-coordinate supports (the crossing/fiber
channels are separate); characteristic-specific asymptotics.
