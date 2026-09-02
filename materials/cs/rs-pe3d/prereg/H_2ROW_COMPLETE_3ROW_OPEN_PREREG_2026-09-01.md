# PREREG-RSPE3D-11-2ROW-COMPLETE-3ROW-OPEN — gate **H-2ROW-COMPLETE-3ROW-OPEN**

Written 2026-09-01 (local date), agent `RsPe3dH2`, after H-DP2-SIZE5 was
frozen, closed FROZEN-CERTIFIED, independently verified by Main, and committed.
This file is path-scoped committed before any parameter-dependent compute.
Binding order: preregistration → path-scoped commit → `campaign.py init` →
byte-identical prereg copy with commit/SHA-256 → controls/censuses → analytic
records → `freeze` → exactly one `close --verdict`. Write only under
`cs/rs-pe3d/`; never modify a frozen directory or root ledger.

The candidate capstone and higher-row program were supplied by Main as derived
and externally checked, not as authority. The two-row statement must be proved
analytically and independently controlled. The higher-row frontier is promoted
only to the strongest theorem actually proved; census alone never becomes a
closed form.

## 1. P1 — complete the two-2-row theory

For any matrix with column rank at most $R$, prove every circuit has size at
most $R+1$: a circuit of size $m$ has rank $m-1$, hence $m-1\le R$. Apply this
to $A\otimes B$ with two 2-row full-rank factors, whose columns lie in
$F^{2\cdot2}=F^4$:

$$|S|\le5.$$

Equality is attainable by every certified size-five support. A six-set cannot
be a circuit because each of its five-subsets is already dependent in ambient
dimension four.

Combine the frozen theorems into one exhaustive list for two $2\times n$ GRS
factors (spark 3):

1. size 3: exactly T-DGE axis fibers;
2. size 4: exactly profile-$(3,3)$ crossings and profile-$(4,4)$
   all-distinct bilinear/Möbius/cross-ratio supports;
3. size 5: exactly the H-DP2-SIZE5 graph/minor templates;
4. no circuit of any other size.

Restate the degree bound with correct scope. For tied spark-$d$ MDS factors, a
circuit of size greater than $d$ has every row/column degree at most $d-1$,
because $d$ cells on one axis contain a proper size-$d$ fiber circuit. More
generally for MDS sparks $d_A,d_B$, a nonfiber circuit larger than both
minimum fiber sizes has A-vertex degree at most $d_B-1$ and B-vertex degree at
most $d_A-1$. Do **not** assert a universal `<=d-1` bound in the unequal case:
a size-$d_B$ B-fiber with $d_B>d_A=d$ is a circuit and is the explicit
exception.

Prove the general ambient bound

$$\boxed{|S|\le r_A r_B+1}$$

for any two-factor product with factor row ranks at most $r_A,r_B$, and state
that equality is attainable (the two-row size-five families furnish examples;
no claim that every factor pair attains it).

## 2. Exact two-row controls

For $2\times4$ and $2\times5$ Vandermonde factor pairs on points $1..n$ over
both GF(7) and GF(13):

- exhaust every six-subset ($\binom{16}{6}=8008$ and
  $\binom{25}{6}=177100$);
- assert zero six-circuits by direct rank/proper-subset checking;
- assert every six-set contains a dependent five-subset;
- retain known size-five ACCEPT supports and assert rank four/all deletions
  independent, demonstrating equality in the $R+1$ bound.

## 3. P2 — open the $2\times n$ by $3\times n$ frontier

Let

$$A_x=(1,x)^{\mathsf T},\qquad
B_y=(1,y,y^2)^{\mathsf T},$$

with distinct points. Then $(d_A,d_B,d)=(3,4,3)$, and product columns

$$h(x,y)=(1,y,y^2,x,xy,xy^2)^{\mathsf T}$$

lie in dimension six, so circuits can have sizes only 3 through 7.

Prove the already-covered layers:

- size 3: exactly A-factor triples at fixed B index, count
  $n\binom n3$, profile $(3,1)$;
- size 4: exactly B-factor quadruples at fixed A index, count
  $n\binom n4$, profile $(1,4)$. There is no nonfiber: a size-four crossing
  would require $d_A+d_B-2=5$, and the all-distinct rank condition is barred
  by $d_A+d_B>d+3$.

For every nonfiber circuit, reuse dual isolation to prove the projection floor
$|\pi_A S|\ge3$, $|\pi_BS|\ge4$. At sizes 5–7, a proper B-fiber quadruple or
A-fiber triple is forbidden, so the correct degree bounds in the convention
where A indices are rows are

$$\deg_A\le d_B-1=3,
\qquad \deg_B\le d_A-1=2.$$

(The reversed wording is a convention hazard; the numerical bounds must be
attached to the factor supplying the potential fiber.)

### Size 5

The crossing window is exactly
$d_A+d_B-2=5$. Prove that every profile-$(3,4)$ five-circuit is a crossing and
that every crossing is a circuit. Its closed count is

$$\boxed{N_{34}^{(5)}=d_Ad_B\binom n3\binom n4
=12\binom n3\binom n4.}$$

Classify every additional measured profile by a finite bipartite
component/degree template and an explicit rank-four/all-four-subsets-independent
predicate. Promote another field-free channel only if an analytic
parameterization, injectivity proof, and exact formula are supplied.

### Sizes 6 and 7

A size-six circuit must have rank five and every five-subset independent. A
size-seven circuit must have rank six and every six-subset independent (all
seven columns are automatically dependent in ambient dimension six). For each
measured profile/template, record these exact minor conditions. Identify
field-free channels and derive closed counts only where a support
parameterization is proved injective; otherwise label the exact determinant
sum OPEN rather than fitting the two-prime counts.

At size seven, explicitly confirm existence/nonexistence and the ambient-bound
edge: a seven-set is necessarily dependent, but it is a circuit only when all
six-deletions are independent.

## 4. Higher-row exhaustive matrix

Use true nested row-major Kronecker columns and exact GF($p$) arithmetic.
Census both:

1. $A=2\times4$, $B=3\times4$ on points 1..4;
2. $A=2\times5$, $B=3\times5$ on points 1..5;

at $p\in\{7,13\}$. Exhaust sizes 3–7 (for the analytic frontier tables,
sizes 5–7 are load-bearing). For every size/profile record:

- total supports swept and exact circuit support set;
- profile counts;
- a label-invariant component/degree template signature (component edge and
  vertex counts, cycle rank, sorted A/B degree sequences);
- number stable across the two primes and number prime-specific as support
  sets;
- direct rank/deletion-minor acceptance.

Assert the size-3/4 formulas and the size-5 crossing support set exactly. For
all sizes 5–7 assert every circuit satisfies projection/degree bounds. Where a
closed construction is proved, generate it independently and require full set
equality, not count equality.

## 5. Controls in both directions

1. **Known-true two-row size-five ACCEPT** showing $R+1=5$ attained.
2. **Two-row six-set REJECT:** every selected six-set has a dependent
   five-subset and cannot be a circuit.
3. **Known-true unequal crossing ACCEPT:** construct one
   $(d_A,d_B)=(3,4)$ crossing of size five and verify exact rank/minimality.
4. **Inserted proper fiber REJECT:** add cells around a known A triple or B
   quadruple so the larger label set contains that dependent proper subset.
5. **Ambient-dimension plant:** choose seven distinct unequal-factor product
   columns, assert rank at most six/dependence automatically, and separately
   test whether deletion minors make it a circuit; dependence alone must never
   be treated as minimality.
6. **Duplicate-factor-column REJECT:** clone one factor column, assert spark
   drops to 2 and a dependent pair appears.
7. **Non-tensor label REJECT:** corrupt one actual product column while
   preserving a label-valid known support; actual rank/minor checks must reject
   even if label logic accepts.
8. **True 2-by-3 Kronecker/concat guard:** product dimension six with exact
   multiplicative coordinates; block concatenation has dimension five and
   must reject.

## 6. Promotion and lifecycle

FROZEN-CERTIFIED for the whole gate only if the two-row capstone is proved and
controlled **and** the requested higher-row sizes 5–7 are analytically complete
by profile/template with every field-free formula and every determinantal
predicate justified. If the higher-row census exposes channels without a
proved closed parameterization/template exhaustion, close
FROZEN-INCONCLUSIVE. The final record must still isolate the fully proved
2-row theorem, state exactly which unequal-row pieces are certified, and give
the strongest exhaustive range/counts and first remaining gap.

All mathematical decisions use Python integers modulo $p$, exact
Gauss–Jordan ranks/determinants, stdlib only. Set
`sys.dont_write_bytecode=True`; no frozen imports; execute with
`PYTHONDONTWRITEBYTECODE=1`. Run one process at nice 15 with
`OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=VECLIB_MAXIMUM_THREADS=1`.
Hard budget: CPU 600 seconds, wall 900 seconds. Any mismatch writes and
preserves `broken_results.json`; disclose/fix the cause and rerun affected
families. No parameter-dependent compute before this prereg commit. If amended
after init, provenance carries both hashes.
