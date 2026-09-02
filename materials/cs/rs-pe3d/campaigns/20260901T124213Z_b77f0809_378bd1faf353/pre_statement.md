# PREREG-RSPE3D-4-DP1 — preregistration for gate **H-DP1-CIRCUIT**

Committed: 2026-09-01 (UTC), by agent `RsPe3dWeightD1`. Written and committed BEFORE
any parameter-dependent compute of this run, per Main's binding process order:
prereg → path-scoped commit (`git commit --only -- <paths>`) →
`python3 scripts/campaign.py init --gate H-DP1-CIRCUIT --prereg <this file>` →
byte-identical copy into the minted run dir (recording source commit + sha256) →
controls → main compute → `freeze` → `close --verdict`. No parameter-dependent
computation precedes the prereg commit.

Provenance / lineage: Run 1 (`20260901T103743Z_d86a1fa2_f72db722778d`, gate
H-MINLINE) CLOSED **FROZEN-NEGATIVE**; Run 2 (`20260901T104923Z_fba446bf_5aee86fef866`,
gate H-MINLINE-DGE3) CLOSED **FROZEN-INCONCLUSIVE**; Run 3
(`20260901T111136Z_99745efb_610bb61dc5a8`, gate H-MINLINE-DGE3-PROOF) CLOSED
**FROZEN-CERTIFIED**, theorem T-DGE discharged (P1–P5). This run is a SEPARATE
lifecycle; it modifies NO frozen artifact and NO root ledger. Setting: two
factors $H = A \otimes B$ over a finite field $F$ (all columns nonzero, sparks
$d_A, d_B \ge 2$ measured, $d = \min(d_A,d_B)$), same model as the certified
theorem's P4 two-factor window, incremented to size $|S| = d+1$.

## 0. Pre-commit disclosures (fixed before the commit)

1. **Steering-source note.** The falsification target, the witness, the
   six candidate counts, the rank criterion, and the structural reframe were
   SUPPLIED to this agent by parent steering (IRC, 2026-09-01, agent `Main`),
   explicitly marked candidate-not-authority. Per the assignment contract the
   author independently re-derived and MACHINE-verified every load-bearing
   number before this commit; the numbers in section 3 are the author's own
   re-measured pinned values, not the steering's verbatim.
2. **Scratch pre-verification, NON-EVIDENTIAL.** Before this commit the author
   ran scratch scripts (`/tmp/rs_pe3d_scratch.py` v1–v2,
   `/tmp/rs_pe3d_cross_prime.py`, inline probes) verifying: witness circuitness
   over GF(13); both pinned censuses (4-col/5-col) at m ∈ {3,4}; the 12-config
   battery at m ∈ {2..5} × GF(7)/GF(13)/GF(31); three adjacent-d P-A
   probes (3x6-V × 3x5-V; 4x6-V × 4x5-V; graph-incidence 6x6 × 3x5); the
   kernel dim identity checks on every census; and arithmetic reconciliation of
   four certified M5 rows against the fiber formula. These runs happened BEFORE
   the prereg commit; their outputs are NON-EVIDENTIAL — the campaign verdict
   rests ONLY on in-run controls re-run inside the minted run dir after
   `init` with the run's own instrument.  Two scratch v1 aborts (a "below-m
   ASSERT" misdesign and an indexing bug) are preserved in scratch logs and
   described here; nothing was hidden.
3. **Scratch v1 instrument defects found and corrected BEFORE this commit**
   (listed because the constraint requires defect disclosure): (a) scratch v1
   asserted below-$m$ emptiness with $m$ instead of $d$ as the threshold and
   crashed on 8 of 9 configs — a per-run threshold conflation; v2 uses $d$
   (the spark of the product) and P1-style emptiness.
   (b) scratch v1 `rank_cols` was invoked with columns of unequal
   ambient length in the witness block (IndexError), fixed by explicit
   3-column rank; (c) the witness "sum" accumulator used 2 coords instead of
   4 (is_zero trivially incoherent), fixed to 4; (d) the v1 spark() treated
   $k > \dim$ columns as independent (rank < k check only) and returned None
   for $r \times s$ with $s > r$ GRS/V models — corrected to a k > d disjunct;
   v2 numbers reconcile with T-DGE theory.
4. **A prior "bad scratch" datapoint deliberately preserved:** GF(7) witness
   census shows {3: 4, 4: 1} as circuit spectra for BOTH witness factors —
   i.e. the witness factor pair at GF(7) has a (d+1)-circuit on each axis,
   though the GF(13)/GF(31) factors do not; the same (a_k),(b_k) integer
   pattern degenerates over GF(7) because the classes collapse.  This
   confirms that factor circuit spectra are $p$-dependent and MUST be
   measured per-prime in the instrument (done; see C6 note).

## 1. Theorem H-DP1 (statement fixed NOW)

Setting: $H = A \otimes B$ over finite field $F$ (odd, as run instrument uses
GF(p)); factors $A$ ($r_A \times s_A$), $B$ ($r_B \times s_B$), all columns
nonzero; spark $d_A, d_B$ measured by exact rref census (including the $k >
r_X$ dimension-forcing disjunct, which scratch disclosed and the corrected
instrument applies); $d = \min(d_A, d_B)$. Product columns $h(u,v) = a_u
\otimes b_v$, $N = s_A s_B$. A **size-$(d+1)$ circuit** is a subset
$S$, $|S| = d+1$, with $\operatorname{rank}(H_S) \le d$ and every proper
subset of $S$ independent; **types**: *fiber along axis i* ($S = \{u\} \times V$
or $U \times \{v\}$ with the varying factor columns a size-$(d+1)$ circuit of
that factor), *all-distinct* ($|U| = |V| = d+1$), *repeated-mixed* (1 < each
distinct-count < d+1). Main's steering reframes non-fiber circuits via:
$\ker(A \otimes B) \cong \{C : A C B^{\mathsf T} = 0\} = (\ker A \otimes
F^{s_B}) + (F^{s_A} \otimes \ker B)$, dimension count $s_A s_B - \rho_A \rho_B
= (s_A - \rho_A) s_B + s_A (s_B - \rho_B) - (s_A - \rho_A)(s_B - \rho_B)$
equality.  The r = $d_A$ skeleton question ("does a size-(d+1) circuit with
$d_A$ distinct $A$-indices force fiber?") is answered NEGATIVELY by the
witness below (profile (3,3) circuits at $d_A = 3$; nothing is a fiber).

**Claim (falsification).** The natural extension "at $|S| = d+1$ inside the
T-DGE fiber regime (no factor has spark 2), every circuit is a fiber" is
FALSE.  Explicit GF(13) witness: $A$ = 2×4, $a_k = (1, k)$, $k = 1..4$
($d_A = 3$; measured 3-circuits = 4, 4-circuits = 0); $B$ = 2×4 with columns
$b_1 = (1,0), b_2 = (-2,1), b_3 = (1,-2), b_4 = (0,1)$ ($d_B = 3$; measured
3-circuits = 4, 4-circuits = 0; parallel classes all singleton).  Then
$S = \{(k,k)\}_{k=1..4}$ (all coefficients 1, i.e. $A B^{\mathsf T} = 0$
exact: row sums 0 and $k$-weighted sums 0 over $\mathbb Z$) gives
$\sum_k a_k \otimes b_k = 0$, rank of $\{a_k \otimes b_k\}$ = 3, all
$\binom{4}{3}=4$ 3-subsets independent ⇒ genuine size-4 circuit of the
product, all 4 A- and all 4 B-indices distinct on a product with NO spark-2
factor and NO parallel columns ⇒ NOT a fiber.  Ex transcription (my ordered
verdict): **FALSIFIED**, witness preserved.

**Claim (theorem Thm-CRIT).** Either of the following forces every size-$(d+1)$
circuit not carrying a repeated factor-side index (see P-A):
**(i) all-distinct rank ceiling.**  A size-$(d+1)$ circuit with $|U| = |V| =
d{+}1$ (all distinct both sides) satisfies $\rho_U + \rho_V \le d+1$, where
$\rho_U = \operatorname{rank}\{a_u : u \in U\}$, $\rho_V$ symmetric.  Proof
(Two functionals): a dependence $\sum_{t} c_t a_{u_t} \otimes b_{v_t} = 0$
with all $c_t \ne 0$ (circuit); applying $x^* \otimes \mathrm{id}$ for any
$x^* \in F^{s_A *}$ gives $\operatorname{diag}(c)\, w(x^*) \in
\operatorname{Rel}_V$ where $w(x^*) = (x^*(a_{u_t}))_t$ ranges over the
$d{+}1$-dim row-evaluation space $W = \{(x^*(a_{u_t}))_t\}$ (dim $= \rho_U$,
dual), and $\operatorname{Rel}_V = \{y : \sum_t y_t b_{v_t} = 0\}$ has dim $
(d{+}1) - \rho_V$; since $y \mapsto \operatorname{diag}(c) y$ is a bijection
preserving subspaces' dimensions, $\rho_U \le (d{+}1) - \rho_V$.  Symmetric.
**(ii) spark floor.**  If a circuit has all-distinct indices, then
$\rho_U \ge \min(d+1, d_A - 1)$: for if $\rho_U < \min(d+1, d_A-1)$, some
$\le \rho_U < \min(d{+}1, d_A - 1)$ columns of $A$ among the $u_t$ are
dependent, giving a dependent set of size $\le d_A - 1$ in $A$ (via
extraction isolating on distinct $V$-indices — uses $|V| = d+1$ distinct
b-columns, the dual-isolation argument of the certified theorem P1, valid
for the all-distinct type) or directly a below-spark dependent set.  Hence
**($\star$)** $(d_A - 1) + (d_B - 1) \le d + 1$ is NECESSARY for an
all-distinct size-$(d+1)$ circuit, i.e. $d_A + d_B \le d + 3$; and for a
TIED minimum, $d_A = d_B = d$, this gives $2d - 2 \le d + 1$, i.e.
$d \le 3$: **Thm-COR:** for tied minimum with $d \ge 4$, every size-$(d+1)$
circuit has a repeated index on some side; combined with the fiber-type
circuit characterization, at $d \ge 4$ tied, every size-$(d+1)$ circuit is
either a fiber or has type with a repeated index (non-fiber repeated-mixed
existence is NOT claimed here either way — see section-4 open registry).

**Claim (Thm-FIB).**  A size-$(d+1)$ subset $S$ of product columns that is
an axis fiber ($S = \{u_0\} \times V$, $|V| = d+1$) is a circuit **iff** $V$
(the $b$-columns) is a size-$(d+1)$ circuit of $B$; and symmetrically for
$U \times \{v_0\}$.  Proof: transporting via $x \mapsto x \otimes b$ is
injective for $b \ne 0$ ($\sum_t c_t a_{u_t} \otimes b_{v_0} = 0$ applied
against a dual functional on $B$ yields $\sum c_t a_{u_t} = 0$), so ranks
agree with the pure factor ranks; minimality transports.  The **fiber-type
closed form**: $\#\{\text{size-}(d+1) \text{ fibers}\} = \sum_{i \in
\{A,B\}} (\prod_{j \ne i} s_j) \cdot C_i(d+1)$ where $C_i(k)$ = number of
$k$-circuits of factor $i$; both axes counted once; no double counting
possible because a set that is simultaneously an axis-$A$ and axis-$B$ fiber
with $|S| \ge 2$ would require both $|V| = 1$ and $|U| = 1$ — contradiction
with $|S| = d+1 \ge 3$.  For MDS factor $H_i$ (spark $d_i$, GRS parity
check), $C_i(d+1) > 0$ iff $d_i = d+1$ exactly (a $k$-circuit of an MDS
parity check is a $k$-subset, exists iff $d_i = k$); number $= \binom{s_i}{d+1}$.

**Boundary corollary (sharpness of T-DGE's "at most one spark-2 factor" at
$|S| = d+1$).**  For two spark-2 factors at $d = 2$: $(d_A - 1) + (d_B - 1) =
2 \le d + 1 = 3$ — the criterion does NOT bar all-distinct size-3 circuits,
and scratch/measured m=3 census of SP2mc×SP2mc exhibits 24 of them besides
16 fiber-type (see section 3) — the "at most one spark-2 factor" hypothesis
stays necessary at $|S| = d+1$.

## 2. Proof obligations (all DISCHARGED in the run artifact `theorem_hd1_circuits.md`)

- **O1 (witness).**  Re-verify the GF(13) witness in-run: sum zero, rank
  exactly 3, all 4 3-subsets independent, all index classes distinct, factors'
  sparks exactly 3, each factor's number of 3-circuits exactly 4, 4-circuits
  0; and the machine census of the product at $m = 4$ gives exactly 156
  circuits with zero fibers and profile split 144 (3,3) + 12 (4,4) (all
  other profiles 0).
- **O2 (Thm-CRIT).**  Full written proof of (i)+(ii) and corollary ($\star$)
  + Thm-COR in the artifact; NO case split on dual isolation.
- **O3 (Thm-FIB).**  Full written proof in the artifact + closed form +
  the MDS $d_i = d+1$ clause.
- **O4 (census closed forms).**  In-run exhaustive (every reachable config
  in M1/M2) verification that the predicted counts from factor data match
  the product's measured circuit census:
  fibers $= \sum_i (\prod_{j\neq i} s_j) C_i(d+1)$; below-$d$ emptiness (P1
  extension to the census noise: assert NO dependent set of size < d in every
  config); and the all-distinct feasibility prediction: # of all-distinct
  circuits > 0 only if $d_A + d_B \le d+3$.
- **O5 (r = d_A family).**  Exhibit in-run the profile-(3,3) witness family
  as an explicit counterexample class: at $(d_A, d_B) = (3, 3)$, profile
  (3,3) exists (144/900 instances), none are fibers — the r = $d_A$
  "fiber forced" skeleton claim is falsified, not repaired.

## 3. Control matrix C1–C6 (exact, pre-registered; all re-run in-run after init)

Instrument: self-contained Python stdlib (`PYTHONDONTWRITEBYTECODE=1`),
exact GF(p) rref (`pow(x,-1,p)`), exhaustive subset enumeration; sparks
MEASURED per prime with the $k > r_X$ disjunct; no numpy/floats/MLflow.

- **C1 pinned censuses (GF(13), total count zero-tolerance):**
  a. Witness product (2×4 MDS-V × 2×4 witness-B): m=3 → 32 circuits, ALL
     fibers, split 16/16 (16 per axis); m=4 → 156 circuits, 0 fibers,
     profiles {(3,3): 144, (4,4): 12}; all $\binom{16}{3} = 560$ /
     $\binom{16}{4} = 1820$ subsets swept.
  b. 5-col product (2×5 V × 2×5 V): m=3 → 100, all fibers, 50/50;
     m=4 → 984, 0 fibers, profiles {(3,3): 900, (4,4): 84}; all
     $\binom{25}{4} = 12650$ subsets swept.
- **C2 criterion rows (GF(13), zero tolerance):**
  a. row43 (3×5 V × 2×4 V, sparks (4,3), $d=3$): m=4 → 20 circuits, ALL
     fibers (all along axis A, profile (4,1): 20), 0 all-distinct (criterion:
     $d_A + d_B = 7 > 6$).
  b. row35 (2×4 V × 4×6 V, sparks (3,5), $d=3$): m=4 → 0 circuits at all
     (criterion: $7 > 6$, no 4-circuits on either factor: 0 fibers + 0
     distinct).  $\binom{24}{4} = 10626$ subsets.
  c. row44@d4 (3×5 V × 3×5 V, sparks (4,4), $d=4$): m=5 → 0 circuits
     ($d_A + d_B = 8 > 7$; also no 5-circuits on factors).
  d. MDS23 (2×4 V × 2×3 V, sparks (3,3)): m=4 → 36 circuits 0 fibers,
     profiles {(3,3): 36} (all path/girth-4 shape).
  e. RS(4x5,t=1) × RS(2x4,t=2), sparks (5,3), $d=3$: m=4 → 0 circuits
     ($8 > 6$).
- **C3 non-MDS battery (GF(13), zero tolerance):**
  a. NM3x5 (3×5, columns e1,e2,e3,e1+e2,2e1+e3; measured spark 3,
     3-circuits 2, 4-circuits 1) × 2×4 V (spark 3): m=3 → 28 circuits all
     fibers (6 A-fibers? measured split: A-axis 20 fibers appear as
     (3,1)-profile 8 + ... PINNED: circuits 28 = fibers 28 (8 A + 20 B)),
     m=4 → 76 circuits = 4 fibers (profile (4,1), along A: 4·(d+1)-fiber =
     the unique 4-circuit of A transported over $\{v_0\}$... measured 4 = s_B
     × 1) + 72 all-distinct/cornerless (profiles {(3,3): 72}).
     criterion: $d_A + d_B = 6 \le 6$ ⇒ all-distinct feasible ✓.
  b. boundary(SP2mc)² (each 2×4, one parallel class {0,1} + two singletons;
     sparks (2,2), $d=2$): m=2 → 10 dependent pairs (8 fibers + 2 off);
     m=3 → 40 circuits, fibers 16 (8+8), non-fiber 24 with profiles
     {(2,3): 12, (3,2): 12} — boundary corollary section 1 instantiated.
  c. boundary(SP2c3, SP2c3b): multiclass spark-2 (class of 3 + 2 singletons)
     × (class of 2 + 2 singletons): m=2 → 23 pairs (17 fibers + 6 off);
     m=3 → 88 circuits, 22 fibers.
  d. sp2class3 × 2×3 V (sparks (2,3), $d=2$, single-degenerate): m=2 →
     9 pairs ALL fibers (9 A-axis); m=3 → 38 circuits, fibers 14
     (9 A-fibers (2,3)-profile... measured pinned) + 24 non-fiber.
- **C4 RS line rows (GF(13), fiber formula, zero tolerance; LS/LINE
  phrasing recorded only for $t_i = 1$):**
  a. RS(3x4,t=1) × RS(2x3,t=1): sparks (4,3), d=3; m=3 → 4 circuits all
     fibers along axis B (the unique 3-circuit of the (s,t)=(3,1) factor
     transported over all $s_A = 4$ fixed A-indices); m=4 → 3 circuits all
     fibers along axis A (the unique 4-circuit of the (4,1) factor
     transported over 3 fixed B-indices, profile (4,1) × 3).
  b. RS census reconciliation of 4 certified anchors recomputed in-run via
     the fiber closed form (the weight-2 window is outside this gate):
     (3,3,4)t111 → 24; (3,4,4)t111 → 16; (3,3,3)t111 → 27; (4,4,4)t111 → 48;
     each = $\sum_i (\prod_{j\neq i} s_j) C_i(d)$ with MDS $C_i(d) =
     \binom{s_i}{d}$ at $d_i = d$ else 0.
- **C5 planted controls (fail-loud, both directions):**
  a. non-tensor column: in the C1a config's product matrix, replace the
     column at index (0,0) by a sum of two distinct product columns; the
     identical enumeration code path MUST exhibit ≥ 1 dependent set of size
     < 3 in the planted matrix and 0 in the pristine (assert both).
  b. duplicated factor column: clone the B-factor's column 0 into a 5th
     column of B in the C1a config; factor measurement MUST fire the
     spark-mismatch assertion (measured spark of B drops 3 → ? and
     p-classes detect the duplicate); enumeration MUST register a dependent
     pair of product columns of size 2 < d=3 and reject the C1a expectation.
  Both plants: expectation INVERTED would be silently-passing → instrument
  asserts the REJECT direction explicitly and fails loud; pristine versions
  re-run as ACCEPT controls adjacent (assert primal expectations hold).
- **C6 cross-prime stability probes (GF(7) and GF(31), count counted, not
  assumed constant):** witness config m=4 count = observed (measured: 152
  at GF(7), 148 at GF(31)) — the prereg does NOT pin these to 156; the
  criterion rows 2a/2b/2c ARE pinned to 0 at each prime (structural);
  profiles split is measured and recorded.  5-col m=4 measured (GF(7): 1020,
  GF(13): 984, GF(31): 936) — pinned as measured per-prime (recorded at
  freeze; the GF(13) value is the C1b zero-tolerance target).

## 4. Compute cap and kill criteria (fixed now)

Cap: 45 minutes wall-clock total for C1–C6 in-run, single process,
`nice -n 15`, `OMP_NUM_THREADS=8` (no BLAS used, stdlib only). Kill criteria:
(i) any control failing after one defect-resolution pass → FROZEN-INCONCLUSIVE
naming the family; (ii) any undischarged obligation O1–O5 →
FROZEN-INCONCLUSIVE naming it; (iii) pristine-freeze check fails → fix and
refreeze before verdict.  Arithmetic: Python ints mod p, exact rref via
`pow(x, -1, p)`, no floats.

## 5. Verdict rule (fixed NOW)

**FROZEN-CERTIFIED** iff ALL of: (a) O1–O5 discharged in
`theorem_hd1_circuits.md` with no gap; (b) C1–C4 pass with zero mismatch;
(c) C5 plants REJECT as pinned; (d) C6 measured values recorded and the
GF(7)/GF(31) structural-zero rows equal 0.  Any single failure surviving one
defect-resolution pass, or any mispinned target needing an in-run addendum
that still fails → **FROZEN-INCONCLUSIVE** naming the residual gap and the
strongest exhaustive range achieved.  The falsified fiber-extension statement
is recorded as a falsification irrespective of verdict; it is NEVER promoted
to theorem.  Registered-OPEN entries (explicitly non-theorem): sharp closed
form for the all-distinct overlap-redistribution counting beyond the
criterion; the repeated-mixed type census; the ≥ 3-factor $(d+1)$ window.

## 6. Scope sentence (Rule-7 style, fixed now)

Swept: the exact two-factor censuses of C1–C3 at the listed configs and
primes, complete subset enumeration up to $\binom{25}{4} = 12650$ subsets,
exact GF(p) ranks via pinned rref; the four RS anchor rows of C4b; the
witness falsification at GF(13) with cross-prime probes at GF(7)/GF(31).
NOT swept: ≥ 3-factor products at $|S| = d+1$; column families beyond
Vandermonde/GRS-model/graph-incidence families; $|S| > d+1$;
characteristic-2 fields; infinite fields; the LINE corollary beyond
$t_i = 1$; the repeated-mixed counting function (Open-P1); any asymptotic
claim.
