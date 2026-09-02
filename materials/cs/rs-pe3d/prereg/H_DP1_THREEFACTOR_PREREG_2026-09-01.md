# PREREG-RSPE3D-9-DP1-THREEFACTOR — gate **H-DP1-THREEFACTOR**

Written 2026-09-01 (local date), agent `RsPe3dH2`, after
H-ALLDISTINCT-PGLCOUNT was frozen, closed FROZEN-CERTIFIED, independently
verified by Main, and committed. This file is path-scoped committed before any
parameter-dependent compute for this run. Binding order: preregistration →
`git commit --only` → `campaign.py init --gate H-DP1-THREEFACTOR` →
byte-identical prereg copy in the minted run directory with source commit and
SHA-256 → exact controls and censuses → proof record → `freeze` →
`close --verdict` with exactly one verdict. No frozen campaign directory is
modified.

The candidate theorem, proof skeleton, and census anchors below were supplied
by Main as derived and externally checked, not as authority. The candidate
contained one explicitly honest completeness gap; this campaign must close it
analytically before promotion and reproduce all anchors independently.

## 1. Fixed theorem and notation

Let

$$H=H_1\otimes\cdots\otimes H_n,
\qquad d_i=\operatorname{spark}(H_i)\ge3,
\qquad d=\min_i d_i,$$

with every factor column nonzero. A support $S$ is a size-$(d+1)$ circuit when
$|S|=d+1$, the corresponding product columns are dependent, and every proper
subset is independent. Coordinate $i$ **varies** on $S$ if its projection has
more than one index.

The target theorem is:

> Every size-$(d+1)$ circuit varies in at most two coordinates. A one-coordinate
> support is an axis fiber over a size-$(d+1)$ factor circuit. A two-coordinate
> support is a nonfiber size-$(d+1)$ circuit of the corresponding two-factor
> product with every other coordinate fixed. There is no genuinely
> three-dimensional size-$(d+1)$ circuit.

For correct non-overlapping count semantics, define:

- $C_i(d+1)$: number of size-$(d+1)$ circuits of factor $H_i$;
- $N^{\rm nf}_{ij}(d+1)$: number of size-$(d+1)$ circuits of
  $H_i\otimes H_j$ whose two coordinate projections both have size greater
  than one (the **nonfiber** two-coordinate population).

The exact count to prove is

$$\boxed{N_n(d+1)=
\sum_i\left(\prod_{k\ne i}s_k\right)C_i(d+1)
+\sum_{i<j}\left(\prod_{k\notin\{i,j\}}s_k\right)
N^{\rm nf}_{ij}(d+1).}$$

The superscript `nf` is load-bearing: if $N_2$ meant all two-factor circuits,
pair fibers would be counted once for every irrelevant companion coordinate
and then again by the first sum. The displayed disjoint formula is the precise
version of the candidate count.

## 2. Proof obligations

### P1 — close the two-factor completeness gap

Let $A\otimes B$ have nonzero columns and sparks $d_A,d_B\ge3$, put
$d=\min(d_A,d_B)$, and let $S$ be a size-$(d+1)$ circuit with full-support
relation coefficients $\gamma_e\ne0$.

1. **Projection floor by dual isolation.** Write $p=|U|$, $q=|V|$. Group the
   relation by its $A$ indices:
   $$\sum_{u\in U}a_u\otimes w_u=0,
   \qquad w_u=\sum_{v:(u,v)\in S}\gamma_{uv}b_v.$$
   If $p<d_A$, the distinct $a_u$ are independent and dual isolation forces
   every $w_u=0$. Every occupied row would then carry a nonzero $B$-relation
   and have degree at least $d_B$. With two rows this costs
   $2d_B>d+1$ cells. Thus either $p=1$ and $S$ is a fiber, or $p\ge d_A$.
   Symmetrically, every nonfiber circuit has $q\ge d_B$.
2. **Rank-sum inequality for repeated as well as distinct profiles.** Enumerate
   the $m=d+1$ cells as edges and form edge-column matrices $A_E,B_E$, including
   repetitions. The relation is
   $$A_E\operatorname{diag}(\gamma)B_E^{\mathsf T}=0.$$
   Since the diagonal is invertible, the row space of
   $A_E\operatorname{diag}(\gamma)$ has dimension $\rho_A$ and lies in the
   relation space of $B_E$, of dimension $m-\rho_B$. Therefore
   $$\rho_A+\rho_B\le d+1.$$
   The projection floor and factor sparks give
   $\rho_A\ge d_A-1$, $\rho_B\ge d_B-1$, hence
   $$d_A+d_B\le d+3.$$
   Since both sparks are at least $d\ge3$, a nonfiber circuit can survive only
   at $d=3$ with $d_A=d_B=3$.
3. **Exhaust the four residual profiles.** At $d_A=d_B=3$ and $|S|=4$, the
   only possible nonfiber profiles are $(3,3),(3,4),(4,3),(4,4)$.
   The frozen H-MIX-H2-CLOSURE theorem identifies every $(3,3)$ circuit as a
   crossing. The $(4,4)$ case is exactly the all-distinct two-factor channel
   and is retained inside $N^{\rm nf}_{ij}$ for arbitrary factors; no GRS
   structure or Möbius normal form is assumed in the general theorem. The
   frozen H-ALLDISTINCT-CROSSRATIO/PGLCOUNT theorems provide that sharper
   bilinear/Möbius classification only when both factors are 2-row GRS
   factors, as in T1–T4. Exclude
   $(3,4)$ analytically: the three $A$ vertices have degrees $(2,1,1)$; the
   rank-sum inequality forces their columns to have rank two and hence a
   one-dimensional full-support relation. The kernel of $A_U\otimes I$ is
   then this relation line tensored with the $B$ space. The two singleton-row
   vectors force its common tensor vector to be proportional to two distinct
   $B$ columns, impossible because spark$(B)=3$ makes every pair independent.
   Exclude $(4,3)$ symmetrically.

Conclude the complete two-factor statement:

- fibers are exactly transported size-$(d+1)$ factor circuits;
- every nonfiber size-$(d+1)$ circuit exists only for $d_A=d_B=d=3$ and is
  exactly crossing profile $(3,3)$ or all-distinct profile $(4,4)$.

This obligation must be proved line by line in-run. If any step fails, the
general theorem cannot be certified.

### P2 — induction through a composite factor

Induct on $n$. Regroup $A=H_1$ and
$B=\bigotimes_{j=2}^n H_j$. Frozen T-DGE gives
$\operatorname{spark}(B)=d_B:=\min_{j\ge2}d_j$ and says every size-$d_B$
circuit of $B$ is an axis fiber.

Apply P1 to $A\otimes B$.

1. A circuit with the $B$ index fixed is a transported size-$(d+1)$ circuit
   of $H_1$, so it varies in one coordinate.
2. A circuit with the $H_1$ index fixed is a size-$(d+1)$ circuit of $B$.
   If $d_B=d$, apply the induction hypothesis to $B$; if $d_B=d+1$, T-DGE
   makes it an axis fiber; if $d_B>d+1$, it is below spark and impossible.
3. A crossing can occur only at $d=d_1=d_B=3$. Its size-three $B$ circuit is,
   by T-DGE, an axis fiber in one coordinate $j\ge2$. Thus the global crossing
   varies only in coordinates $1,j$.
4. An all-distinct circuit can occur only at $d=d_1=d_B=3$. Its four involved
   $B$ columns have rank at most two by the rank-sum inequality and
   $\rho_A\ge2$. Hence every three-subset is dependent and, by spark$(B)=3$,
   is a size-three circuit. T-DGE makes each such triple an axis fiber. Two
   triples share two distinct columns; two axis lines sharing two points have
   the same varying coordinate and the same fixed coordinates. Hence all four
   $B$ columns lie on one axis line, so the configuration varies in exactly
   one coordinate of $B$ and at most two globally.

Prove that tensoring by the fixed nonzero columns is injective, so the induced
support on coordinates $1,j$ is the corresponding genuine two-factor circuit,
not merely a projection lookalike.

### P3 — disjoint exact count

Partition every circuit by its unique set of varying coordinates, whose size
is one or two by P2. One varying coordinate $i$ contributes
$(\prod_{k\ne i}s_k)C_i(d+1)$. Two varying coordinates $i,j$ contribute
$(\prod_{k\notin\{i,j\}}s_k)N^{\rm nf}_{ij}(d+1)$. Prove these classes are
pairwise disjoint and exhaust all circuits, yielding the displayed formula.

**Promotion rule.** FROZEN-CERTIFIED only if P1 closes the predecessor's
repeated-mixed gap with no hidden completeness assumption, P2–P3 are complete,
all listed exhaustive support sets equal their constructed families exactly,
every numeric anchor is reproduced, and every planted REJECT fires. Otherwise
close FROZEN-INCONCLUSIVE and name the first exact proof/control gap plus the
strongest exhaustive range.

## 3. Exact control matrix

All mathematical decisions use Python integers modulo $p$, exact Gauss–Jordan
rank, and actual nested row-major Kronecker columns. Every circuit check asserts
full-set dependence and every proper-subset independence. The constructed
family is generated independently from factor circuits and two-factor
**nonfiber** circuit censuses, then compared as a support set to the full
three-factor product census.

### T1 — three $2\times3$ Vandermonde factors over GF(7)

Use columns $(1,x)^\mathsf T$, $x=1,2,3$, in all three factors ($N=27$).
Exhaust all $\binom{27}{3}=2925$ triples and $\binom{27}{4}=17550$
four-subsets. Pin:

- exactly 27 size-three circuits, all fibers, 9 per axis;
- exactly 81 size-four circuits;
- size-four profile split
  $(1,3,3):27$, $(3,1,3):27$, $(3,3,1):27$;
- zero supports with all three coordinates varying;
- pair formula $81=3$ unordered pairs $\times3$ fixed third indices
  $\times9$ two-factor crossing circuits.

Assert exact equality of the measured 81-support set and the constructed
pairwise family.

### T2 — three $2\times4$ Vandermonde factors over GF(7)

Use columns $(1,x)^\mathsf T$, $x=1,2,3,4$ ($N=64$). Exhaust all
$\binom{64}{4}=635376$ four-subsets. Independently census the two-factor
reference and pin exactly 152 nonfiber circuits, split 144 crossing plus 8
all-distinct. Pin the three-factor census:

- 1824 size-four circuits;
- profiles $(1,3,3),(3,1,3),(3,3,1)$: 576 each;
- profiles $(1,4,4),(4,1,4),(4,4,1)$: 32 each;
- zero genuinely three-dimensional supports;
- $1824=3\times4\times152$.

Assert exact measured/constructed support-set equality.

### T3 — a factor carrying size-$(d+1)$ circuits, with unequal sparks

Use over GF(7):

- $M_{34}$: $3\times4$ Vandermonde columns $(1,x,x^2)^\mathsf T$,
  $x=1,2,3,4$, spark 4, with its unique four-subset a size-four circuit;
- two $V_4$: $2\times4$ Vandermonde factors, spark 3.

Exhaust all 635376 four-subsets of $M_{34}\otimes V_4\otimes V_4$. Pin:

- 16 one-axis fibers of profile $(4,1,1)$;
- 576 crossings of profile $(1,3,3)$;
- 32 all-distinct pairs of profile $(1,4,4)$;
- total 624, no other profile and zero genuinely three-dimensional supports;
- the $M_{34}\otimes V_4$ pair has zero nonfiber size-four circuit, while
  $V_4\otimes V_4$ has 152.

Assert exact measured/constructed support-set equality. This is both an
unequal-spark control and a direct test of the separate single-factor summand.

### T4 — unequal sizes/sparks independent cross-check

Use $M_{34}\otimes V_3\otimes V_4$ over GF(7), with $V_3$ points 1..3.
Exhaust all $\binom{48}{4}=194580$ four-subsets. Pin exactly:

- 12 fibers of profile $(4,1,1)$;
- 144 pair circuits of profile $(1,3,3)$, because
  $V_3\otimes V_4$ has 36 nonfiber crossings and the $M_{34}$ index has four
  choices;
- total 156, zero other profiles and zero genuinely three-dimensional
  supports.

Assert exact measured/constructed support-set equality.

### T5 — controls in both directions

1. **Known-true ACCEPT.** Materialize a crossing support in
   $V_3\otimes V_3$ and tensor every column by one fixed nonzero column of a
   third $V_3$; assert rank 3 and every triple independent.
2. **Non-tensor REJECT.** On that same labeled support, replace the first true
   three-factor column by itself plus the column with the same first two
   indices and a distinct third index. Assert actual rank rises to four while
   the label-only pairwise constructor still predicts the support.
3. **Genuine-3D REJECT.** In $V_4^{\otimes3}$ use diagonal cells
   $(t,t,t)$ for all four indices. Assert all three projections vary but the
   actual four columns have rank four, so the support is not a circuit and is
   absent from the constructed family.
4. **Duplicated-column REJECT.** Clone one factor column and assert its spark
   drops from 3 to 2 with a dependent pair.
5. **True 3-factor Kronecker guard.** For three 2-row columns assert nested
   Kronecker dimension 8 and exact multiplicative coordinates. Plant block
   concatenation (dimension 6) and a dropped-third-factor product (dimension
   4); assert both differ and reject.
6. **Formula-semantics REJECT.** Deliberately use all pair circuits instead of
   nonfiber pair circuits in T3. Assert the resulting multiset double-counts
   the $M_{34}$ single-axis fibers and does not equal the disjoint theorem
   family.

## 4. Arithmetic, budget, defects, and kill criteria

- Exact mathematical arithmetic only: Python integers modulo 7; inverses by
  `pow(x,-1,p)`; no floating-point mathematical decisions, NumPy, or external
  algebra. Timing does not enter a theorem decision.
- `sys.dont_write_bytecode=True`; no frozen import; execute with
  `PYTHONDONTWRITEBYTECODE=1`.
- One process at `nice -n 15`, with
  `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=VECLIB_MAXIMUM_THREADS=1`.
- Hard budget: CPU 600 seconds, wall 900 seconds. Abort rather than raise it.
- Every positive claim has a known-true ACCEPT and a planted REJECT in-run.
- Any assertion mismatch is fail-loud and writes `broken_results.json` before
  nonzero exit. Preserve every broken artifact. Disclose any instrument or
  prereg defect, repair the causal issue without weakening pins, and rerun all
  affected families.
- No parameter-dependent compute may occur before this file's path-scoped
  commit. If this prereg is amended after init, provenance records both hashes.
- Do not infer a general theorem from census. The proof obligations P1–P3,
  especially the $(3,4)/(4,3)$ exclusion and the unequal-$d_B$ induction
  branches, are mandatory even if every anchor passes.
