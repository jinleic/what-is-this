# PREREG-RSPE3D-20-C3-34-RESIDUAL-POINT — gate **H-C3-34-RESIDUAL-POINT**

Written 2026-09-04 UTC by agent `RsSparkClose`, immediately after
H-SPARK2-CONVERSE was frozen and closed FROZEN-CERTIFIED (run
`20260904T040116Z_00df70c5_9bf664de2def`). This file is path-scoped committed
before any new parameter-dependent computation. Binding order:
preregistration -> path-scoped commit -> `campaign.py init` -> byte-identical
in-run copy with source commit and SHA-256 -> standalone instrument -> controls
before search -> checkpoints -> analytic/experimental record -> freeze ->
exactly one close verdict. No frozen run and no root `cs/RESULTS.md` is edited.

## 1. Why this is the shortest remaining Rs frontier

The frozen C3 direction concerns the size-$(r_A r_B-1)$ pencil layer. At
$(r_A,r_B)=(3,3)$, two coprime $(2,2)$ equations have intersection length
8, exactly the support size; fixed-component alternatives were classified by
the size-8 residual gate. The first strictly new case is $(3,4)$: the
function space

$$V=H^0(\mathbf P^1\!\times\!\mathbf P^1,\mathcal O(2,3))$$

has dimension 12, a C3 support has 11 points, and two coprime equations have
intersection length 12. Thus a new residual point appears. This gate binds
the universal residual-length statement separately from a finite GF(17)
existence/minimality experiment. A finite-field search can certify only its
registered family; it cannot promote a universal converse or a
characteristic-zero claim.

This gate also carries the existing all-distinct existence caveat: the frozen
law says no all-distinct circuit exists below $K=r_A+r_B$ and classifies the
size-$K$ layer **when present**. It does not assert that a PGL population, a
C3 population, or the product spark equals an all-distinct threshold when the
corresponding incidence set is empty.

## 2. Exact geometric hypotheses and universal lemma

Work over a field $F$ and then over its algebraic closure. Let $S$ be 11
distinct $F$-rational affine points of $\mathbf P^1\times\mathbf P^1$ whose
$x$- and $y$-projections are both injective. Let

$$W_S=\{H\in V:H|_S=0\}.$$

The **pencil hypotheses** are:

1. $\dim_F W_S=2$, equivalently the 12-by-11 evaluation matrix has rank 10;
2. a basis $F_0,G_0$ of $W_S$ has no common irreducible curve component over
   the algebraic closure (basis-independently, the pencil has no fixed
   component);
3. every point of $S$ is a reduced point of the scheme-theoretic intersection
   $Z=V(F_0,G_0)$ (local intersection multiplicity one).

**Lemma R1 (universal, proof obligation rather than experimental inference).**
Under these hypotheses, Bezout on $\mathbf P^1\times\mathbf P^1$ gives

$$\operatorname{length} Z=(2,3)\cdot(2,3)=2\cdot3+3\cdot2=12.$$

The residual scheme $Q=\operatorname{Res}_S(Z)$ therefore has length exactly
one. Since $Z$ and the 11 reduced rational points are defined over $F$, $Q$ is
a single reduced $F$-rational projective point (possibly outside the affine
evaluation grid). This is only a residual-length/rationality lemma.

The rank-10 hypothesis makes $S$ dependent with one-dimensional column
relation space. Circuit minimality is a separate condition:

$$S\text{ is a circuit}\iff
\operatorname{rank}E_{S\setminus\{s\}}=10\quad\text{for every }s\in S.
\tag{M}
$$

The gate does **not** assume that no fixed component alone implies (M).
Finding one support satisfying the three pencil hypotheses but failing (M)
is a valid negative outcome for that sufficiency conjecture.

## 3. Registered GF(17) construction family

The experiment uses exact arithmetic over GF(17) and constructs two
bidegree-$(2,3)$ equations in factored form

$$F=L_M Q_R,\qquad G=L_N Q_T.$$

Here $L_M=(cx+d)y-(ax+b)$ is the irreducible $(1,1)$ graph of a PGL map
($ad-bc\ne0$), and

$$Q_R=x(d_0+d_1y+d_2y^2)-(n_0+n_1y+n_2y^2)$$

is the irreducible $(1,2)$ graph of a degree-two rational map from the
$y$-line: numerator and denominator are projectively coprime and their
maximum degree is exactly two. All coefficient tuples are canonically scaled
by their first nonzero entry. A candidate pair requires $M\ne N$ and
$R\ne T$; unique factorization and the unequal factor classes then prove
that $F,G$ share no component.

The instrument deterministically builds 96 valid canonical PGL factors and
96 valid canonical degree-two factors from the 31-bit LCG

$$s_{k+1}=1103515245s_k+12345\pmod {2^{31}},$$

with pool seed `0xC33417`. It then inspects the first **1,000,000 distinct
unordered compatible equation pairs** generated with search seed
`0x34C317`, skipping repeats, identical equations, and shared-factor pairs.
This cap, seeds, pool sizes, prime, and ordering are fixed; no adaptive
extension is allowed.

For each pair, the affine common-zero set is computed by two exact routes on
the full 17-by-17 grid: direct products `F==G==0`, and the union of the four
factor-pair intersections. The routes must agree whenever a set is promoted.
A **qualifying 12-set** $Z_0$ must have exactly 12 distinct affine common
zeros and 12 distinct $x$-values and 12 distinct $y$-values. Because the
factors are coprime and Bezout's total intersection length is 12, these 12
points exhaust the projective intersection and are reduced; there is no
unseen point or multiplicity.

For every qualifying $Z_0$, all twelve 11-point deletions
$S_q=Z_0\setminus\{q\}$ are checked. Promotion requires, separately for each
$S_q$:

1. rank exactly 10 **before** any deletion/minimality check (so
   $W_{S_q}=\langle F,G\rangle$ really is a pencil);
2. one-dimensional, full-support relation space;
3. every one of its eleven 10-point deletions has rank exactly 10.

Thus each promoted base set receives 12 rank-before-deletion checks and 132
10-deletion checks. Rank-deficient 11-sets and rank-10 but nonminimal 11-sets
are retained as distinct negative classes with the first exact witness. A
counterexample is never suppressed merely because another deletion is a
circuit.

## 4. Registered outcomes and claim boundary

- **R1** is universal and must be proved algebraically as stated; it does not
  depend on finding a GF(17) example.
- **E-positive:** at least one qualifying GF(17) 12-set exists and all of its
  checked 11-deletions are classified exactly. The record may certify only
  the completed one-million-pair registered experiment, including its exact
  population and whether every qualifying deletion was a circuit.
- **E-negative:** a qualifying rank-10 11-set fails a 10-deletion rank. This
  falsifies no-fixed-component sufficiency on the registered GF(17) family;
  preserve the support, equations, residual point, ranks, and failing
  deletion.
- **E-empty:** no qualifying 12-set occurs in the fixed family. Report exact
  emptiness as a registered-family result, but use FROZEN-INCONCLUSIVE for the
  existence/minimality frontier. Do not infer global nonexistence.

No finite outcome proves a universal C3 converse, an extension-field claim,
or a characteristic-zero claim. A universal converse could be promoted only
from a separately written proof that discharges (M), never from sample size.

## 5. Controls before the search

All controls use the same exact rank/evaluation/factor code as the search.

1. Dual rank implementations agree on a 12-column tensor-product Vandermonde
   basis and on every later promoted support.
2. **Circuit ACCEPT:** ten independent ambient basis vectors plus their
   full-support sum form an 11-vector rank-10 circuit; every deletion has rank
   10.
3. **Full-rank promotion REJECT:** eleven independent vectors have all
   deletions of rank 10 but rank 11 before deletion, so a deletion-only bug
   must not promote them.
4. **Nonminimal dependent REJECT:** an 11-vector set containing a planted
   duplicate subcircuit is dependent but fails deletion minimality.
5. **Geometric ACCEPT:** the hard-coded frozen `(3,3)` complete-intersection
   support `((1,5),(2,6),(3,1),(4,4),(5,2),(6,3),(7,8),(8,7))` over GF(13)
   is checked standalone (no frozen import) as rank 7 with every deletion rank
   7.
6. Valid PGL and coprime degree-two factors are accepted; a singular PGL tuple
   and a numerator/denominator common-factor tuple are rejected.
7. Distinct factor pairs are accepted as no-common-component by exact
   canonical factor identities; planted shared-`L` and shared-`Q` pairs are
   rejected.
8. Direct-grid and four-factor-intersection common-zero routes agree on a
   planted pair; corrupting one factor evaluation makes the equality fail.
9. A 12-point set with a repeated coordinate is rejected by the all-distinct
   gate.

The results file and checkpoint are written immediately after these controls,
before the first search pair, so a genuinely launched partial run has an
exact Stage-F record.

## 6. Runtime, preservation, and verdict rule

`run_c3_34_residual.py` is standalone and imports no frozen campaign code.
Exact GF(p) integers only; no floats, NumPy, or SymPy. Set
`sys.dont_write_bytecode=True` before other imports and launch with
`PYTHONDONTWRITEBYTECODE=1`; all five thread caps equal 1; `nice -n 10` and
`RLIMIT_CPU=(5400,5400)` are asserted. Hard caps are 5400 CPU seconds and
6000 wall seconds. Checkpoint every 10,000 distinct pairs and on every
qualifying 12-set. Any failure writes a uniquely named broken-results file;
no file or failed attempt is deleted. A defect is logged and its affected
stage rerun from the beginning without weakening assertions.

A terminal **FROZEN-CERTIFIED** verdict requires R1 proved and a nonempty
completed registered experiment with every candidate classified under the
fixed rules. **FROZEN-NEGATIVE** is allowed if the completed experiment
contains a qualifying rank-10 nonminimal 11-set, with its full witness, while
R1 remains a certified subresult. **FROZEN-INCONCLUSIVE** applies to an empty,
incomplete, or failed registered experiment. Exactly one verdict may be
closed; nothing inferred from a finite search is phrased universally.
