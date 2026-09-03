# PREREG-RSPE3D-15-GRS-GENERAL-SIZE — gate **H-GRS-GENERAL-SIZE**

Written 2026-09-02 UTC by agent `RsPe3dH2`, after H-33ROW-SIZE8-RESIDUAL was
frozen, closed FROZEN-CERTIFIED (run
`20260902T041030Z_64e0c4ef_b5de16e49f8a`), independently verified by Main, and
committed. This file is path-scoped committed before any new
parameter-dependent computation. Binding order: preregistration -> commit ->
`campaign.py init` -> byte-identical in-run copy with source commit and
SHA-256 -> controls before censuses -> analytic record -> `freeze` -> exactly
one `close --verdict`. Write only under `cs/rs-pe3d/`; do not edit any frozen
run or root ledger.

The candidate claims below were derived analytically before registration and
cross-checked against the owner's derivation, but they are not authority. The
run must re-verify every proof obligation line by line and confirm every
positive claim in exact arithmetic, each with a known-true ACCEPT and a
planted REJECT. Counts outside proved constructions remain census data.

## 1. Setting, all-distinct stratum, fallback predicates

Let the two factors be GRS/Vandermonde matrices of sizes $r_A\times n_A$ and
$r_B\times n_B$ over $F_p$,

$$a_x=(1,x,\ldots,x^{r_A-1})^{\mathsf T},\qquad
b_y=(1,y,\ldots,y^{r_B-1})^{\mathsf T},$$

with distinct evaluation points inside each factor and arbitrary nonzero
column multipliers. Nonzero multipliers rescale product columns and change no
rank, circuit, or count below, so they are suppressed. The product column is
the row-major Kronecker evaluation

$$h(x,y)=a_x\otimes b_y=(x^i y^j)_{0\le i<r_A,\ 0\le j<r_B}\in F^{r_Ar_B}.$$

For a support $S$ of grid cells write
$E_S=\bigl[h(x,y)^{\mathsf T}\bigr]_{(x,y)\in S}$ (the $r_Ar_B\times|S|$
evaluation matrix). **All-distinct** means coordinate-injective: the
$x$-coordinates of $S$ are pairwise distinct and the $y$-coordinates are
pairwise distinct. Crossings, shared-fiber, and other tied profiles are
separate strata and are never folded into the all-distinct layer.

The exact fallback predicate (C) is

$$S\text{ circuit}\iff \operatorname{rank}E_S=|S|-1\text{ and every
one-column deletion of }E_S\text{ has rank }|S|-1.$$

The ambient cap is $|S|\le r_Ar_B+1$ (proved in gate 7 without the GRS
hypothesis); at $|S|=r_Ar_B+1$ dependence is automatic and (C) is exactly the
nonvanishing of all $|S|$ deletion minors; at $|S|=r_Ar_B$ minimality needs no
deletion beyond rank corank one plus every deletion independent.

## 2. Theorem T — transport, with proof obligations

**T.** For all $r_A,r_B\ge2$ and all-distinct $S$: $S$ is a circuit of
$A\otimes B$ if and only if $\operatorname{rank}E_S=|S|-1$ and every proper
subset of $S$ is independent. Moreover the same transport holds for the
tied strata through the obligations below; ranks of arbitrary supports are
computed on actual Kronecker columns throughout.

Proof obligations to be verified line by line in-run:

(a) **Column form.** The column of the true product matrix at $(x,y)$ is a
nonzero scalar multiple of $h(x,y)$; multipliers never affect rank or
minimality.

(b) **Kronecker-vs-slice identity.** $(A\otimes B)\operatorname{vec}(\Gamma)=0
\iff A\,\Gamma\,B^{\mathsf T}=0$, with $\operatorname{vec}$ the same row-major
order as $h$; row-major Kronecker coordinates and dimensions asserted against
a block-concatenation impostor.

(c) **Minimality transport.** Independence of a sub-support and the deletion
predicate (C) transport between $A\otimes B$ and $E_S$; the ambient cap and
the corank-one reading of (C) at sizes $r_Ar_B$ and $r_Ar_B+1$ follow.

T is verified numerically on every support the run tests; property (b) is
additionally stressed by the concat guard and a corrupted-column plant.

## 3. C1 — the unified all-distinct law at size $r_A+r_B$ (Möbius/PGL)

**C1 (forward direction — proved).** Let $S$ be all-distinct with
$|S|=r_A+r_B\ge5$. If the pairs of $S$ lie on one nondegenerate Möbius graph
— there are $(c_3,c_2,c_1,c_0)\ne0$ with $\Delta=c_0c_3-c_1c_2\ne0$ and
$y=-\dfrac{c_2x+c_0}{c_3x+c_1}$ on $S$, finite poles excluded by distinctness
inside $Y$ — then $S$ is dependent of the expected rank
$r_A+r_B-1$. This forward implication is proved in general for arbitrary
finite $X,Y$ and primes.

**C1 per-prime equality (registered; verified numerically in-run).** At every
registered anchor and prime the run must verify by exhaustive support-set
comparison:

- **(eq-a)** the set of all-distinct dependent $k$-sets, $k=r_A+r_B$, equals
  the set of Möbius-graph $k$-sets; and
- **(eq-b)** its cardinality equals the PGL sum

$$N=\sum_{M\in\mathrm{PGL}(2,p)}\binom{|D_M|}{k},\qquad
D_M=\{x\in X: M(x)\in Y\},$$

enumerating $M$ by the canonical normalization first nonzero of
$(c_3,c_2,c_1,c_0)$ to one, $|\mathrm{PGL}(2,p)|=p(p^2-1)$; graphs emitted
twice must be asserted against map redundancy ($p-1$ scalars share a class).
These equalities are registered per prime (§6), not as a universal closed
form. The run must also confirm zero all-distinct dependent sets at sizes
below $r_A+r_B$ at the anchors; those zeros are census anchors (C1-pre).

**C1 converse — SCOPED.** Whether every all-distinct dependent
$(r_A+r_B)$-set must lie on a single Möbius graph for arbitrary $X,Y,p$ is an
open lead; (eq-a) at the anchors is numerical support, not a proof. An
elliptic $(2,2)$ construction with $h^0$ residual $2r_A+2r_B-4$ is consistent
with the converse (no counterexample at large factors), but no general
converse is claimed, and no verdict of this gate depends on it.

**C1 emptiness below $r_A+r_B$ — SCOPED beyond the anchors.** Lead: a slice
nullspace device $\lambda_j=w_jq(x_j)$ with $w_j=\prod_{l\ne j}(x_j-x_l)^{-1}$
and $\deg q\le k-r_A-1$ generalizes the certified $r_A=3$ device; a
double-sided $y$-argument forcing $R'\circ R=\mathrm{id}$ on $k$ points closes
only when $k>\deg R\cdot\deg R'+1$ under $r_Ar_B-3r_A-3r_B+7\le0$ — a
condition the run checks at every registered anchor. This boundary is
recorded; claims beyond the registered anchors are not made.

## 4. C2 — the top two sizes (proved, pure linear algebra)

**C2 (i), size $r_Ar_B$.** An $S$ of size $r_Ar_B$ is a circuit iff
$\operatorname{rank}E_S=r_Ar_B-1$ and every one-column deletion has full rank
$r_Ar_B-1$. Proof: corank one gives one independence; deletions full rank
exclude every proper dependence; converse immediate. Equivalently: $E_S$ has
a one-dimensional nullspace and all $r_Ar_B$ deletion minors are nonzero.

**C2 (ii), size $r_Ar_B+1$.** Every $S$ of size $r_Ar_B+1$ is dependent
(ambient cap). $S$ is a circuit iff every $(r_Ar_B)$-subset of $S$ is
independent, equivalently all $r_Ar_B+1$ deletion minors are nonzero; a
single singular deletion exhibits a proper dependence and fails minimality.
This specializes gate 9's predicates (13)/(14):
- size-9 at $(r_A,r_B)=(3,3)$: $r_Ar_B=9$; predicate (C2)(i) is exactly
  registered predicate (13) with its $(2,2)$-curve determinant reading;
- size-10 at $(3,3)$: $r_Ar_B+1=10$; (C2)(ii) is exactly predicate (14) (all
  ten $9\times9$ deletion determinants nonzero).

Gate 7's size-5 five-minor predicate and the tied three-row size-9/10
predicates are thereby re-expressed as C2 instances; the run must confirm the
re-expression by recomputing the certified census numbers at the registered
grids.

## 5. C3 — the pencil layer at size $r_Ar_B-1$ (one direction + bootstrap)

**C3 (proved direction).** Let $|S|=r_Ar_B-1$, $\dim W_S=2$ in the sense that
the vanishing ideal of the evaluation grid restricted to bidegrees
$(<r_A,<r_B)$ has pencil-like nullspace spanned by two independent equations,
and the pencil has no fixed component. Then the base locus of the pencil is
exactly $S$, which must hence be a circuit of $A\otimes B$. At
$(r_A,r_B)=(3,3)$: two independent bidegree-$(2,2)$ equations without a
fixed component meet in $(2,2)\cdot(2,2)=8=|S|+1$ points, leaving exactly
one residual grid point per pencil; the certified gate-10 dichotomy applies
and the size-8 support is a genuine circuit. Concretely the pipeline is: fit
the exact pencil, compute its common-component content via the gate-10
engine (Sylvester/GCD over $F_p[x]$ and $F_p[y]$), assert no common factor
of positive degree in either variable, verify base locus $=S$ pointwise,
and conclude circuitness by (C).

**C3 bootstrap.** Confirm the $(3,3)$ size-8 complete-intersection law
against gate 10's certified counts: 560 complete intersections over GF(11)
and 416 over GF(13) in the $n=8$ matching sweep, with the rank-7 pencil
counts 2192/1344 and the common-$(1,1)$ noncircuit counts 1632/928
recomputed by the same instruments where swept.

**C3 converse — SCOPED.** The converse fails at $(3,3)$: gate 10 certified
shared-$(1,1)$-component pencils whose base locus strictly contains $S$. The
no-fixed-component characterization of the whole slice is therefore a scoped
open lead, not a claim; the direction proved above is verdict-bearing, the
characterization is not.

## 6. Registered anchors (owner instances, exact reproduction required)

Field: GF(13); points $1..n$; factors exact Vandermonde. Additional prime
p=11, then p=7, replicating each anchor as a per-prime equality claim.

1. $(r_A,r_B,n)=(2,4,6)$, $k=r_A+r_B=6$: direct census of all-distinct
   6-subsets: dependent $=2=$ PGL sum; all-distinct 5-subsets dependent $=0$.
2. $(2,4,7)$, $k=6$: all-distinct size-6 circuits $=42=$ PGL sum.
3. $(3,4,7)$, $k=7$: all-distinct 6-subsets dependent $=0$; all-distinct
   size-7 circuits $=2$, both Möbius graphs, PGL sum $2$.
4. $(3,3)$ cross-checks against the frozen gates: gate 9 §5.2 all-distinct
   size-6 circuits at $n=6$: 12/4/2 over GF(7)/GF(11)/GF(13); at $n=7$:
   588/72/42; gate 9 (13)/(14) counts at $4\times4$: 0/16 at sizes 9/10 and
   8/0/16/0/66/0/16 at sizes 4..10 under all three primes; gate 10's 560/416
   CI counts.
5. C2 grids: $4\times4$ sizes 9–10 over GF(7), GF(11), GF(13); $5\times5$
   GF(7) sizes 9–10 with pinned census 36378/181960.

Blunt sweep scale is bounded: $n=7,k=6$ is $35{,}280$ eligible all-distinct
subsets; $k=7$ is $41{,}160$; all budgets below are sized for far larger.

## 7. Controls (must run before any parameter-dependent census)

Every ACCEPT has a paired REJECT in the same process:

1. **Pure-tensor ACCEPT/REJECT.** A known-true circuit (anchor-verified) is
   accepted by (C); a corrupted non-tensor column (one coordinate replaced by
   another column's value under identical labels) must raise rank and be
   rejected — the label/PGL predicate must not accept it.
2. **Duplicated-column spark REJECT.** Clone one factor column; the spark of
   the duplicated matrix drops from $r_A+r_B$ to $r_A+r_B-1$, asserted
   through an explicit dependent pair at the new spark; the uncloned column
   set must stay at spark $r_A+r_B$.
3. **Concat-vs-Kronecker guard.** On an $r_A\times r_B$ factor pair assert
   true row-major Kronecker dimension and coordinates; the block
   concatenation impostor has the wrong dimension and coordinates.
4. **Discriminating ambient plant.** At size $r_Ar_B+1$ assert
   $\operatorname{rank}=r_Ar_B$ exactly (not merely $\le$) and exhibit a
   per-subset witness; a dependent-but-nonminimal set must be rejected by
   (C2)(ii) through one failed deletion. Repeats the gate-10 audit flag.
5. **Non-Möbius plant.** An all-distinct $(r_A+r_B)$-set chosen off every PGL
   graph (verified by the $k_M$ sums and a direct graph test) must be
   independent: $\det E_S\ne0$ in square case or rank $=|S|$; and at $(2,2)$
   the known wrong-coefficient REJECT $q^2(q-2)$ and selected-pole REJECT of
   gate PGL are replayed.
6. **C2 witness plants.** A corank-one size-$r_Ar_B$ set with one singular
   deletion must be REJECTED; a size-$(r_Ar_B+1)$ set with all deletions
   nonsingular must be ACCEPTED as a circuit.
7. **C3 plants.** A no-fixed-component pencil with base locus exactly $S$ is
   ACCEPTED; a shared-$(1,1)$-component pencil (gate-10 type) is REJECTED
   from circuitness while its dependence is confirmed.

## 8. Stage plan for `run_grs_general_size.py`

- **A — transport library:** build $E_S$ as row-major Kronecker of factor
  evaluation matrices; every rank used in a claim is cross-checked by two
  independent rank algorithms (batch RREF and incremental-basis); obligation
  (b) identity spot-checked on random sparse $\Gamma$.
- **B — PGL layer:** canonical enumeration (normalize first nonzero of
  $(c_3,c_2,c_1,c_0)$ to 1), assert $p(p^2-1)$ maps, per-map $D_M$, $k_M$,
  sums $\binom{k_M}{k}$.
- **C — anchor sweeps (Stage F first):** per prime, per
  $(r_A,r_B,n)$: all-distinct dependent counts and circuit counts by actual
  rank plus deletion witnesses at anchor sizes; reproduce §6 exactly,
  primes in order 13, 11, 7.
- **D — C2 instances:** the grids of §6 item 5; every claimed circuit
  confirmed by deletion-witness minors, every claimed noncircuit by one
  exhibited failure.
- **E — C3 pencils:** $(3,3)$ size-8 pencil classification at the registered
  instances; cross-validate against gate 10's 560/416, 2192/1344, 1632/928.
- **F — control battery:** §7 in full before any of C/D/E censuses; per-stage
  checkpoints; any defect freezes the affected family, preserves broken
  output, and reruns that family from its beginning.

## 9. Arithmetic, runtime discipline, budgets, defects

Exact GF(p) integer arithmetic only: Python ints mod $p$, inverses by
`pow(x,-1,p)`; no floats, no NumPy/sympy, no external algebra. `sys.dont_write_bytecode=True`
first, before any import; process launched with `PYTHONDONTWRITEBYTECODE=1`
and `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=VECLIB_MAXIMUM_THREADS=NUMEXPR_NUM_THREADS=1`.

Launch inside the run directory: `nice -n 10 env PYTHONDONTWRITEBYTECODE=1
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 python3
run_grs_general_size.py` (Main's standing order caps nice at 10; the
instrument asserts `os.nice(0)==10`, overriding the template's 15). Before
the first long stage the process sets
`resource.setrlimit(resource.RLIMIT_CPU, (5400, 5400))` soft then hard.
Budget caps: **5400 CPU seconds / 6000 wall seconds**, checked by every
exhaustive family via the shared `check_budget`; no adaptive extension.

Any defect is disclosed in `defect_log.json`, broken output preserved under a
distinct name, and every affected family rerun from its beginning. No
assertion weakened, no count silently repinned. Kill and close per policy on
any failed control, anchor mismatch, enumeration inconsistency, or budget
violation.

## 10. Promotion rule and verdict

**FROZEN-CERTIFIED** iff: obligations (a)–(c) verified; C1 forward and the
PGL-sum equality hold at every registered prime; C2 (i)/(ii) proved and
confirmed on all registered grids including the gate-9 re-expressions; C3's
proved direction holds and the $(3,3)$ bootstrap matches gate 10's certified
counts; every anchor in §6 reproduced exactly; and the full control battery
passes with every planted REJECT firing. Otherwise **FROZEN-INCONCLUSIVE**,
naming the first failed proof or control while certifying inside the record
every part that did pass (gate-9 precedent). Scoped gaps — C1 converse, C1
general emptiness, C3 no-fixed-component characterization — never block the
verdict; they are recorded as open leads exactly as scoped in §§3–5.

Certified on promotion: Theorem T's transport with obligations (a)–(c); C1's
forward Möbius-graph law and its per-prime PGL-sum equalities; C2 at sizes
$r_Ar_B$ and $r_Ar_B+1$ with the gate-9 predicate re-expressions; C3's proved
direction plus the $(3,3)$ bootstrap certified counts. Not certified: C1's
converse; all-distinct emptiness below $r_A+r_B$ beyond the registered
anchors; a no-fixed-component characterization of the C3 layer; arbitrary
$X,Y$ closed forms beyond the PGL sum; sparks not registered here; three or
more factors; extension fields.
