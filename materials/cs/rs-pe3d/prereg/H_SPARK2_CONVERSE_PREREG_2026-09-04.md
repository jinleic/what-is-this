# PREREG-RSPE3D-19-SPARK2-CONVERSE — gate **H-SPARK2-CONVERSE**

Written 2026-09-04 UTC by agent `RsPe3dConverse`, immediately after
H-GRS-M3-GENERAL was frozen, closed FROZEN-CERTIFIED (run
`20260904T025704Z_7c5edb95_267bbf870b1d`). This file is path-scoped committed
before any new parameter-dependent computation. Binding order:
preregistration -> commit -> `campaign.py init` -> byte-identical in-run copy
with source commit and SHA-256 -> controls before censuses -> analytic record
-> `freeze` -> exactly one `close --verdict`. Write only under `cs/rs-pe3d/`;
do not edit any frozen run or root ledger.

## 1. Setting (framework of the frozen Theorem X, gate H-MIX-CROSSING)

$A$: $r_A\times s_A$, $B$: $r_B\times s_B$ over GF($p$), $p$ odd, nonzero
columns; sparks $d_A,d_B\ge2$ measured (minimum dependent-column-set size);
product columns $h(u,v)=a_u\otimes b_v$ indexed by the grid
$[s_A]\times[s_B]$; a support's **profile** is (distinct $A$-indices, distinct
$B$-indices); $m_* = d_A+d_B-2$; $\mathrm{Circ}_X(k)$ = $k$-circuits of
factor $X$; $C_X(k)$ = their count. A **spark-2 factor** has
$d_X=2\iff$ it possesses a proportional column pair; such pairs are exactly
its 2-circuits, so $C_X(2)$ = number of proportional pairs. The frozen
Theorem X classified size-$m_*$ profile-$(d_A,d_B)$ circuits as exactly the
crossings $S(R,i_0,J,j_0)$ when $d_A,d_B\ge3$, leaving the branches
$d_A=2$ and/or $d_B=2$ (permutation-shaped possibilities) explicitly open.
This gate closes them.

## 2. Candidate classification (independently derived; to be verified)

Write a two-row support on the proportional pair $\{u,u'\}$ ($a_{u'}=\lambda
a_u$, $\lambda\ne0$) and a $d_B$-subset $J$ of $B$-columns as
$S(T)=\{u\}\times T\ \cup\ \{u'\}\times(J\setminus T)$ with
$\emptyset\ne T\ne J$.

**Theorem S2 (converse, both spark-2 branches).** Let $d_A=2$, $d_B\ge2$,
$m_*=d_B$. The size-$m_*$ profile-$(2,d_B)$ circuits of $H=A\otimes B$ are
EXACTLY the supports $S(T)$ above with $J\in\mathrm{Circ}_B(d_B)$ and
$\emptyset\ne T\ne J$. Each is a genuine circuit with 1-dimensional relation
space $\gamma_{u,v}=\tau\,c_J(v)$ on $T$, $\gamma_{u',v}=-\tau\lambda^{-1}
c_J(v)$ on $J\setminus T$ (where $c_J$ is the $B$-circuit relation), and the
count is
$$N_{(2,d_B)} = C_A(2)\,C_B(d_B)\,\bigl(2^{d_B}-2\bigr).$$
Symmetrically for $(d_A,2)$. For $(2,2)$ this gives
$C_A(2)C_B(2)(2^2-2)=2C_A(2)C_B(2)$, matching the frozen count.

**Proof sketch (full proof in the record).** Profile $(2,d_B)$ with $m_*=d_B$
cells forces a two-row support with each $B$-index used exactly once. Row
relation: $A\Gamma B^{\mathsf T}=a_u(B\alpha)^{\mathsf T}+
a_{u'}(B\beta)^{\mathsf T}$; if $a_u,a_{u'}$ independent this forces
$B\alpha=B\beta=0$, impossible at sizes $<2d_B$ unless a row is empty
(profile drops). If $a_{u'}=\lambda a_u$: $B(\lambda\alpha+\beta)=0$, so the
combined support is a dependent $B$-set of size $d_B=d_B$, i.e. a
$B$-circuit $J$ with $\lambda\alpha+\beta\propto c_J$, giving exactly the
family above. Relation space is 1-dimensional with full support, so every
proper subset is independent. □

**Dichotomy (registered alternative outcomes).** The crossings for fixed
(pair $\{u,u'\}$, $J$) are exactly the subfamily $|T|\in\{1,d_B-1\}$
($2d_B$ supports). Therefore:
- **Crossing-only** holds $\iff$ $2^{d_B}-2=2d_B\iff d_B\in\{2,3\}$: at
  $(2,2)$ and $(2,3)$ (and symmetrically $(3,2)$) every profile-$(2,\cdot)$
  circuit IS a crossing; the converse boundary closes with no extra family.
- **Explicit extra families** exist $\iff d_B\ge4$: exactly
  $2^{d_B}-2-2d_B$ non-crossing permutation-shaped circuits per
  (proportional pair, $B$-circuit), and Theorem X's "exactly crossings"
  clause is FALSE there; the corrected count replaces
  $2d_B\,C_A(2)C_B(d_B)$ (crossings only) by
  $(2^{d_B}-2)\,C_A(2)C_B(d_B)$.

## 3. In-run obligations (exact GF(p), dual ranks, deletion witnesses)

- **(O1)** On every enumerated size-$m_*$ profile-$(2,d_B)$ (and
  $(d_A,2)$, $(2,2)$) support: circuit predicate (rank $m_*-1$, every
  one-deletion rank $m_*-1$) and relation space dimension 1 with the exact
  proportional-to-$c_J$ coefficient pattern of Theorem S2.
- **(O2)** Completeness sweep: ALL size-$m_*$ supports with profile
  $(2,d_B)$ on the registered factors are enumerated (the restricted
  two-row grid makes this exhaustive) and classified; the population must
  equal the formula count exactly, per prime.
- **(O3)** Counterclassification check: every enumerated circuit that is not
  a crossing must be exhibited as a partition family member with
  $2\le|T|\le d_B-2$; every crossing must have $|T|\in\{1,d_B-1\}$; at
  $d_B\in\{2,3\}$ the non-crossing population must be ZERO.
- **(O4)** Independence below $m_*$: all profile-$(2,d_B)$ supports of sizes
  $2..m_*-1$ on the two-row grid are independent (exhaustive at registered
  sizes).

## 4. Registered instances (factors built in-run; $C_A(2)=1$ by construction)

- **(2,4) branch**: $A$ = 3-row GRS on $\{1..6\}$ with column $a_{6}:=a_{1}$
  duplicated (spark 2, one proportional pair, $C_A(2)=1$); $B$ = 3-row GRS on
  $\{1..5\}$ ($d_B=4$, $C_B(4)$ computed in-run). Primes 13, 11, 7.
  Predicted population $C_B(4)\cdot(2^4-2)=14\,C_B(4)$; crossings
  $8\,C_B(4)$; extra $6\,C_B(4)$ per prime.
- **(2,3) branch**: same $A$; $B$ = 2-row GRS on $\{1..5\}$ ($d_B=3$).
  Predicted: crossing-only, population $6\,C_B(3)$; non-crossing count ZERO
  at all primes.
- **(2,2) branch**: $A$ as above; $B$ = 3-row GRS with one duplicated column
  ($d_B=2$, $C_B(2)=1$). Predicted population $2C_A(2)C_B(2)=2$; both
  circuits are crossings (the two diagonals). Primes 13, 11, 7.
- **Symmetric branch (4,2)**: $A$/$B$ roles swapped from the (2,4)
  instance; primes 13, 11.
- **Independence below $m_*$** (O4): exhaustive on the two-row grids of the
  (2,4) and (2,3) instances at all registered primes.

All counts are pinned by the formula of section 2 BEFORE the sweep (the
formula is the registered prediction); a mismatch aborts the gate.

## 5. Control battery (before censuses; same exact GF(p) path)

1. **Crossing ACCEPT**: the frozen crossing construction $S(R,i_0,J,j_0)$
   built on the spark-2 factor ($R$ = the duplicated pair) is accepted by the
   circuit predicate.
2. **Extra-family ACCEPT** (counterclassification positive): a partition
   support with $2\le|T|\le d_B-2$ at the (2,4) instance is accepted as a
   circuit with the exact coefficient pattern of Theorem S2.
3. **Planted REJECT (sub-circuit)**: every one-deletion of the ACCEPT
   supports is independent (deletion witnesses recorded).
4. **Planted REJECT (wrong profile)**: a row-fiber support
   (profile $(1,d_B)$) of size $d_B$ is DEPENDENT but must be rejected from
   the profile-$(2,d_B)$ classification by its profile.
5. **Planted REJECT (non-proportional pair)**: a two-row support on an
   independent column pair is independent (rank $m_*$).
6. **Corrupted-column REJECT**: replacing one factor column value in the
   ACCEPT support's relation raises the rank and breaks circuitness.
7. **Formula cross-check**: in-run computed $C_A(2)$, $C_B(d_B)$ by
   exhaustive factor-subset census, multiplied per the formula, compared to
   the sweep population (O2) — two independent routes.

## 6. Stage plan for `run_spark2_converse.py`

- **A** library (factor builders incl. duplicated-column GRS, dual ranks,
  nullspace, circuit predicate, crossing constructor, partition enumerator).
- **F** control battery (before censuses).
- **B** factor spectra: sparks, $C_A(2)$, $C_B(d_B)$ by exhaustive census.
- **C** branches (2,4), (2,3), (2,2), (4,2): exhaustive sweeps, (O1)-(O4),
  formula pins, primes in order 13, 11, 7 (13, 11 for the symmetric
  instance).
- **G** wrap up.

## 7. Arithmetic, runtime discipline, budgets, defects

Exact GF(p) integer arithmetic only; no floats, no NumPy/sympy;
sys.dont_write_bytecode=True first; PYTHONDONTWRITEBYTECODE=1; all five
thread caps = 1; nice -n 10 asserted; RLIMIT_CPU soft = hard = 5400 before
the first long stage. Hard cap **5400 CPU s / 6000 wall s**, checked per
family; no adaptive extension. Defects disclosed in `defect_log.json`,
broken output preserved under a distinct name, affected family rerun from
its beginning; no assertion weakened; prereg amendments follow the A1/A2
transparency rule (both hashes, never claimed unamended).

## 8. Promotion rule and verdict

**FROZEN-CERTIFIED** iff: Theorem S2 is registered as proved with its
in-run obligations verified (every enumerated support classified, population
= formula at every registered cell, crossing-only verified at $d_B\in\{2,3\}$
with zero non-crossing circuits, extra families exhibited and counted at
$d_B\ge4$), the full control battery passes with every planted REJECT
firing, and the independence-below-$m_*$ sweeps are clean. Otherwise
**FROZEN-INCONCLUSIVE**, naming the first failed obligation while certifying
inside the record every part that did pass. Nothing incomplete is promoted.
