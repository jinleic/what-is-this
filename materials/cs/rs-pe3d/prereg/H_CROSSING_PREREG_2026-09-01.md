# PREREG-RSPE3D-5-CROSSING — preregistration for gate **H-MIX-CROSSING**

Committed: 2026-09-01 (UTC), by agent `RsPe3dCrossing`. Written and committed BEFORE
any parameter-dependent compute of this run, per Main's binding process order:
prereg → path-scoped commit (`git commit --only -- <paths>`) →
`python3 scripts/campaign.py init --gate H-MIX-CROSSING --prereg <this file>` →
byte-identical copy into the minted run dir (recording source commit + sha256) →
controls → main compute → `freeze` → `close --verdict`. No parameter-dependent
computation precedes the prereg commit.

Provenance / lineage: Run 1 (H-MINLINE) FROZEN-NEGATIVE; Run 2 (H-MINLINE-DGE3)
FROZEN-INCONCLUSIVE; Run 3 (H-MINLINE-DGE3-PROOF, theorem T-DGE P1–P5)
FROZEN-CERTIFIED; Run 4 (`20260901T124213Z_b77f0809_378bd1faf353`,
H-DP1-CIRCUIT) FROZEN-CERTIFIED. This is a SEPARATE lifecycle; it modifies NO
frozen artifact and NO root ledger. Setting: two factors $H = A \otimes B$ over
GF($p$), all factor columns nonzero, sparks $d_A, d_B \ge 2$ measured by exact
rref, $d = \min(d_A, d_B)$; subset size window $|S| = d+1$ plus the exactly
determined mixed window $|S| = d_A + d_B - 2$.

## 0. Pre-commit disclosures (fixed before the commit)

1. **Steering-source note.** The candidate "CROSSING (mixed-profile) channel"
   theorem, the obligation list, and every anchor number below were SUPPLIED to
   this agent by parent steering (agent `Main`, 2026-09-01), explicitly marked
   candidate-not-authority. Supplied anchors: witness totals 152/148/156/148 at
   $p \in \{7,11,13,17\}$ with $(3,3)$-profile 144 at each; 2×5 pair total
   984 = 900 mixed + 84 all-distinct, 0 fibers; 2×6 pair mixed 3600 =
   $9\binom{6}{3}^2$; 3×4 MDS-parity pair $m{=}5 \to 0$, $m{=}6 \to 16$ all
   $(4,4)$; concat anchors 24/16/27/48 (prior-certified). The author
   independently re-derived the construction and all proofs from the frozen
   kernel identity $\ker(A \otimes B) = (\ker A \otimes F^{s_B}) + (F^{s_A}
   \otimes \ker B)$; the theorem in §1 is the author's own sharpened statement
   (bijection with forced-circuit converse, the $d_A,d_B \ge 3$ injectivity
   clause, the $(2,2)$ 2-to-1 halving, and the classified $d_A = 2$ XOR
   $d_B = 2$ singleton cases are author additions). **No scratch or
   parameter-dependent compute was run before this commit**; all numbers below
   are pinned targets from steering + prior certified runs and are verified
   in-run only. Nothing was pre-measured to make the pins safe; a failed pin is
   a possible falsification of a supplied anchor and will be reported as such.
2. **Instrument inheritance.** Helpers are copied from the certified run-4
   instrument (`run_hd1.py`): exact GF($p$) rref via `pow(x,-1,p)`, spark with
   the $k > $ ambient-dim disjunct, `Factor` with zero-column assert, census
   with exhaustive below-$d$ emptiness pass and kernel-dimension identity
   assert. Run 3's defect class (block concatenation instead of Kronecker
   product, `vec[off+i] = f.cols[idx][i]`) is guarded by a dedicated control
   (T8).
3. **Known-lattice note.** The mixed count formula $N_{\text{mix}} =
   d_A d_B\, C_A(d_A) C_B(d_B)$ reproduces all four supplied mixed anchors
   (144 = 3·3·4·4; 900 = 3·3·10·10; 3600 = 9·20·20; 16 = 4·4·1·1). The
   all-distinct residuals (12/84, and 8/4/12/4/1?"(4,4)" profiles) have NO
   closed form (run-4 Thm-CRIT only bounds them); they are pinned only as
   Class-A anchors, never promoted to theorem.

## 1. Theorem H-MIX-CROSSING (statement fixed NOW)

**Setting.** $A$ ($r_A \times s_A$), $B$ ($r_B \times s_B$) over GF($p$), $p$
odd, all columns nonzero; sparks $d_A = \operatorname{spark}(A)$, $d_B =
\operatorname{spark}(B)$ measured; $d = \min(d_A,d_B)$; product columns $h(u,v)
= a_u \otimes b_v$, $N = s_A s_B$; profile of a column set $S$ = $(\#\text{
distinct }A\text{-indices}, \#\text{distinct }B\text{-indices}) = (n_A, n_B)$.
Factor circuit counts $C_X(k)$ = number of $k$-circuits of factor $X$ (an
exact GF($p$) measurement). Kernel identity (frozen, run-3/4):
$\{C : A C B^{\mathsf T} = 0\} = (\ker A \otimes F^{s_B}) + (F^{s_A} \otimes
\ker B)$ with dimension identity $s_A s_B - \rho_A\rho_B = (s_A - \rho_A)s_B +
s_A(s_B - \rho_B) - (s_A-\rho_A)(s_B-\rho_B)$.

**Construction.** For a $d_A$-circuit $R$ of $A$, $i_0 \in R$, a $d_B$-circuit
$J$ of $B$, $j_0 \in J$: let $c$ = the (unique up to scale) circuit relation of
$A$ on $R$ (full support), $\delta$ = that of $B$ on $J$; scale $\delta$ so
that $c_{i_0} + \delta_{j_0} = 0$ (possible: $c_{i_0} \ne 0$, $\delta_{j_0}
\ne 0$). Then $C = c\,e_{j_0}^{\mathsf T} + e_{i_0}\delta^{\mathsf T}$
satisfies $A C B^{\mathsf T} = 0$ and has exact support
$S(R,i_0,J,j_0) := (R \setminus \{i_0\}) \times \{j_0\} \;\cup\; \{i_0\}
\times (J \setminus \{j_0\})$: the **crossing cell** $(i_0,j_0)$ is where both
summands land and cancel.

**Theorem X (crossing channel, mixed-profile circuits).** Let $d_A, d_B \ge 2$
and $m_* = d_A + d_B - 2$. The circuits of $H = A \otimes B$ with
$|S| = m_*$ and profile exactly $(d_A, d_B)$ are precisely the supports
$S(R,i_0,J,j_0)$ with $R \in \mathrm{Circ}_A(d_A)$, $i_0 \in R$, $J \in
\mathrm{Circ}_B(d_B)$, $j_0 \in J$; each carries the explicit full-support
relation of the Construction, has a 1-dimensional relation space (rank $m_*-
1$: the span is $U \otimes \langle b_{j_0}\rangle + \langle a_{i_0}\rangle
\otimes W$ with $U = \operatorname{span}(R \setminus\{i_0\})$, $W =
\operatorname{span}(J\setminus\{j_0\})$, both containing the opposite extreme
vector, so the intersection is exactly $\langle a_{i_0} \otimes b_{j_0}\rangle$),
and the parameterization is:
  - **injective (count $= d_A d_B\, C_A(d_A) C_B(d_B)$)** when $d_A, d_B \ge
    3$ (heavy row/column identification: the $i_0$-row carries $d_B - 1 \ge 2$
    cells, all other rows 1; symmetric for $j_0$);
  - **injective** when exactly one of $d_A, d_B$ equals 2 (the light side's
    singleton cell identifies $j_0$, resp. $i_0$; formula as written);
  - **exactly 2-to-1** when $d_A = d_B = 2$: $(R,i_0,J,j_0) \equiv
    (R,i_1,J,j_1)$, so the count is $2\,C_A(2) C_B(2) = $ formula$/2$.
**Converse (forced circuits).** Conversely every circuit of $H$ of size
$d_A + d_B - 2$ with profile $(d_A,d_B)$ arises this way: for any
decomposition $C = C_1 + C_2$ with $k_1$ nonzero $\ker A$-columns and $k_2$
nonzero $\ker B$-rows, the support bound $|\mathrm{supp}\,C| \ge k_1 d_A +
k_2 d_B - k_1 k_2$ (spark floors, disjoint columns/rows, $k_1k_2$ overlap cap)
plus the profile $(d_A,d_B)$ force $(k_1,k_2) = (1,1)$; the profile forces
$\mathrm{supp}(c) \ni i_0$, $|\mathrm{supp}(c)| = d_A$ (else profile breaks or
spark floor breaks), symmetrically $j_0 \in \mathrm{supp}(\delta)$; size
$d_A+d_B-2$ forces the cancelled center; and circuitness of $S$ plus $c \in
\ker A$ with $|\mathrm{supp}(c)| = d_A$, $\mathrm{supp}(c) \ni i_0$ forces
$R = \mathrm{supp}(c)$ to be a $d_A$-circuit (a proper dependent subset would
let a kernel combination shrink the support below the spark floor), and
symmetrically $J$.
**Field independence.** The count depends on the field only through $C_A(d_A)$
and $C_B(d_B)$: it is an exact polynomial identity in the factor circuit
spectra, valid simultaneously over every GF($p$); the all-distinct population
at the same window has no such closed form and is $p$-dependent (measured:
totals $152/148/156/148$ with $(3,3) = 144$ fixed).
**Reach / coherence with run-4 Thm-CRIT.** The crossing channel produces
$|S| = d+1$ circuits **iff** $d_A + d_B = d + 3$, exactly the tightness of the
run-4 necessary condition $d_A + d_B \le d+3$; the mixed channel is the
explicit saturation of that bound. When $d_A + d_B > d + 3$ the channel lands
at $m_* = d+2$ or beyond (verified: $3{\times}4$ MDS pair lands at 6 = $d+2$
with $m{=}5$ exactly 0).
**Remarks (recorded, not counted).** (a) *L-family*: for $i_0 \in R$, $j_0
\notin J$ with $b_{j_0} \notin \operatorname{span}(J)$, $S = (R \times
\{j_0\}) \cup (\{i_0\} \times J)$ is a circuit of size $d_A + d_B - 1$ with
profile $(d_A, d_B + 1)$; out of every window swept here; counting it is
OPEN. Non-cancelled centers ($i_0 \in R$, $j_0 \in J$, $c_{i_0} +
\delta_{j_0} \ne 0$) are never circuits (they contain the crossing circuit).
(b) *Regrouping corollary* (≥3 factors): $H = H_1 \otimes H_2 \otimes \cdots
\otimes H_n = A \otimes B_{\mathrm{comp}}$ with $B_{\mathrm{comp}} = H_2
\otimes \cdots \otimes H_n$; by certified P1 (T-DGE) $\operatorname{spark}
(B_{\mathrm{comp}}) = \min_i d_i$, and Theorem X applies with $C_B(d_B)$ =
composite-B circuit counts — measured only, no individual-factor closed form.
Certified scope of this gate: the two-factor statement; the regrouping
corollary is verified on one composite instance, not certified in general.
(c) Fiber circuits (profile $(d_A,1)$, $(1,d_B)$; run-4 Thm-FIB) and
all-distinct circuits (run-4 gate) are separate populations, not affected.
(d) At $d = 2$ boundaries the statement sharpens as in the 2-to-1 clause; the
spark-2 degenerate world itself stays out of scope.

## 2. Proof obligations (all to be DISCHARGED in `theorem_crossing.md`)

- **O1 (construction + size law + circuitness).** Written proof that
  $S(R,i_0,J,j_0)$ is a circuit of size $d_A + d_B - 2$, profile $(d_A,d_B)$,
  1-dim relation space; PLUS in-run explicit numeric instance at GF(13)
  ($R = \{0,1,2\}$, $i_0 = 0$; $J = \{0,1,2\}$, $j_0 = 1$ on the witness
  factors): relations in-kernel, scaled center zero, full-support sum zero on
  the 4 product columns, rank $= 3 = m_* - 1$, all 3-subsets independent.
- **O2 (bijection + closed form).** Proof of the parameterization count with
  the three injectivity clauses (incl. 2-to-1 at $(2,2)$) and of the converse;
  in-run **set equality** (not just cardinality) between the census's
  profile-$(d_A,d_B)$ population and the constructed $\{S(R,i_0,J,j_0)\}$ at
  GF(13) for: witness (144), 2×5 (900), 2×6 (3600), 3×4 pair (16),
  2×3V×witnessB (36), composite demo ($9X$).
- **O3 (field independence).** In-run: witness factor spectra $C_A(3) =
  C_B(3) = 4$, 4-circuits 0, at $p \in \{7,11,13,17,31\}$; mixed $(3,3) = 144$
  at each; all-distinct residual $(4,4)$ varies (recorded 8/4/12/4/p31
  measured) — pinned totals Class-A. Written proof that the mixed count is a
  function of factor spectra only.
- **O4 (reach / coherence).** Witness tightness $3+3 = 6 = d+3 \Rightarrow$
  lands at $|S| = 4 = d+1$ (census) ; $3{\times}4$ pair $8 = d+4 > d+3
  \Rightarrow$ 0 circuits at $m = 5$ (exhaustive $\binom{16}{5} = 4368$) and
  16 at $m = 6$ (exhaustive $\binom{16}{6} = 8008$, all $(4,4)$, set equality);
  written proof of the iff.
- **O5 (multi-column exclusion).** Written proof of the support bound $k_1
  d_A + k_2 d_B - k_1 k_2$ and its case analysis forceing $(k_1,k_2) = (1,1)$
  at the counted size; corroborated in-run by census profile purity (no mixed
  profile other than $(d_A,d_B)$ at $m_*$; none at $d+1$ beyond
  $(d,d)$/all-distinct/$m_*$-crossing overlap).
- **O6 (converse).** Written proof (decomposition case analysis + forced
  circuits); machine side is O2's set equalities at all six configs.

## 3. Control matrix (all re-run in-run after `init`; exact GF(p), stdlib only)

Class **T** = theorem-pinned, zero-tolerance (failure ⇒ verdict per §5).
Class **A** = anchor-pinned, zero-tolerance (steering-supplied or
prior-certified; failure after one disclosed defect-resolution pass ⇒ §5).
Global asserts on every census: kernel-dim identity; exhaustive below-$d$
emptiness; fibers/fibers-axis counts; profile bookkeeping; relation-space
dimension $= m-1$ on set-checked populations.

- **T1 (ACCEPT block, GF(13)).** Witness: $A$ = 2×4, $a_k = (1,k)$, $k=1..4$;
  $B$ = 2×4 columns $(1,0), (-2,1), (1,-2), (0,1)$: sparks $(3,3)$;
  $\mathrm{circ}_A = \mathrm{circ}_B = \{3\!:4\}$; 4-circuits 0 both; $B$
  parallel classes all singleton; kernel-dim identity; O1 explicit instance
  asserts.
- **T2/A (witness census GF(13)).** $m{=}4$: circuits 156 [A: prior-certified],
  fibers 0 [T], profile $(3,3) = 144$ [T] with **set equality** vs
  constructed [T], profile $(4,4) = 12$ [A], no other profile [T];
  $\binom{16}{4} = 1820 + 120$ below-d swept. $m{=}3$: 32 circuits, all
  fibers 16/16 [A], 680 subsets swept.
- **T3 (crossing-mechanism plants, GF(13)).**
  - **T3a (cancellation-sensitivity REJECT).** For EVERY constructed crossing
    relation (all 144 at witness, all at 2×5, 2×6, 3×4, 2×3V): perturbing any
    one support coefficient destroys the dependence (rank jumps $m_*{-}1 \to
    m_*$); assert pristine rank $m_*-1$ AND relation-space dim 1 AND full
    support for each constructed relation.
  - **T3b (spectrum-corruption quantitative row).** $A''$ = 2×3 V × witness
    $B$, $d=3$: census $m{=}4$ total exactly 36 [T: $3 \cdot 1 \cdot 4 \cdot 3$;
    all-distinct impossible ($s_{A''}=3 < 4$), fibers 0 (no 4-circuits)],
    set equality of the 36.
  - **T3c (duplicated factor column REJECT, prior C5b).** Clone $B$ col 0 into
    a 5th column: factor measurement fires (spark $< 3$ or multi-column
    pclass) AND product census yields $\ge 1$ dependent pair (size $2 < d$).
  - **T3d (non-tensor column REJECT, prior C5a).** In witness product, replace
    column $(0,0)$ by sum of two distinct product columns: pristine has 0
    non-fiber dependent triples at $m{=}3$, planted has $\ge 1$ (assert both
    directions).
- **T4 (2×5 pair GF(13)).** $m{=}4$: total 984 [A: prior-certified], fibers 0
  [T], $(3,3) = 900$ [T] with set equality, $(4,4) = 84$ [A], no other
  profile [T]; $\binom{25}{4} = 12650 + 300$ swept. $m{=}3$: 100, all fibers
  50/50 [A].
- **T5 (2×6 pair GF(13)).** $m{=}4$ full census $\binom{36}{4} = 58905$,
  below-d $\binom{36}{2} + \binom{36}{3} = 630 + 7140$: $(3,3) = 3600$ [T]
  with set equality, fibers 0 [T], no profile other than $(3,3),(4,4)$ [T];
  $(4,4)$ **measured and recorded, not pinned**; total recorded.
- **T6 (3×4 MDS pair GF(13)).** Factors `vand_rows(3,4,1,13)` both: sparks
  $(4,4)$ [T], $C_A(4) = C_B(4) = 1$, 5-circuits 0 [T: MDS parity check:
  $k$-circuits iff $k = 4$]; $m{=}5$: 0 circuits, 0 dependent [T], 4368
  subsets; $m{=}6$: exactly 16 circuits, all profile $(4,4)$, 0 fibers, set
  equality vs $4\cdot4\cdot1\cdot1$ [T]; $C(16,6) = 8008$ subsets; below-d
  $k\in\{2,3\}$: 680 subsets, 0 dependent [T]; no other profile [T].
- **T7 (cross-prime field independence).** Witness config at $p \in
  \{7,11,13,17,31\}$: factor spectra $\{3\!:4\}$, 4-circuits 0 [T]; $m{=}4$
  profile $(3,3) = 144$ [T], fibers 0 [T]; totals $152/148/156/148/148$ [A];
  $(4,4)$ residuals $8/4/12/4/4$ [A-derived, recorded].
- **T8 (concat guard).** True nested-Kronecker 3-factor censuses reproduce
  $(3,3,4)t(1,1,1) \to 24$, $(3,4,4)t(1,1,1) \to 16$, $(3,3,3)t(1,1,1) \to
  27$, $(4,4,4)t(1,1,1) \to 48$ [A: prior-certified]; the run-3 defect path
  (block concatenation, `vec[off+i] = f.cols[idx][i]`, dim $r_1{+}r_2{+}r_3$
  vs true $r_1 r_2 r_3$) is built and asserted structurally different
  (dimension 7 vs 12 at $(3,3,4)$) with its census count recorded.
- **T9 (regrouping demo, GF(13)).** Composite $B_{\mathrm{comp}} = H_2 \otimes
  H_3$ for $(s,t) = ((3,1),(4,1))$ [i.e. $2{\times}3$V $\otimes$ 3×4 MDS]:
  measured spark $= 3$ [T: P1-IH], spectra recorded; census of $A'' \otimes
  B_{\mathrm{comp}}$ at $m{=}4$ over $\binom{36}{4}$ subsets: $(3,3) = 9 \cdot
  C_{B_\mathrm{comp}}(3)$ with set equality [T], fibers $= 3 \cdot
  C_{B_\mathrm{comp}}(4)$ profile $(1,4)$ [T], no other profile [T], total
  asserted to the measured-spectrum prediction [T].

## 4. Compute cap and kill criteria (fixed now)

Cap: **60 minutes wall-clock** single process, `nice -n 15`,
`OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1`, `PYTHONDONTWRITEB​YTECODE=1`,
Python ints mod $p$, exact rref `pow(x,-1,p)`, no floats/numpy. Kill criteria:
(i) any Class-T or Class-A control failing after ONE disclosed
defect-resolution pass (defects disclosed, broken outputs preserved, affected
families re-run) → **FROZEN-INCONCLUSIVE** naming the family; (ii) any
undischarged obligation O1–O6 → FROZEN-INCONCLUSIVE naming it; (iii) cap
exceeded → FROZEN-INCONCLUSIVE with partial outputs preserved; (iv) any pin
mismatch that survives → FROZEN-INCONCLUSIVE naming the anchor (a failed
steering anchor is reported as a falsification of that anchor, never silently
corrected). Heartbeat via `campaign.py heartbeat` if compute runs long.

## 5. Verdict rule (fixed NOW)

**FROZEN-CERTIFIED** iff ALL of: (a) O1–O6 discharged in
`theorem_crossing.md` with no gap; (b) every Class-T control passes exactly;
(c) every Class-A control passes exactly; (d) recorded-but-not-pinned rows
written to `controls_results.json`. Any single failure surviving the one
defect-resolution pass → **FROZEN-INCONCLUSIVE**, naming the residual gap and
the strongest exhaustive range achieved. Certified scope: Theorem X as stated
in §1 (profile $(d_A,d_B)$, size $d_A{+}d_B{-}2$, injectivity clauses,
converse, field independence, reach iff), at the swept configs/primes only.
The all-distinct counts, the L-family counting function, and general ≥3-factor
closed forms are NOT certified by this gate.

## 6. Scope sentence (Rule-7 style, fixed now)

Swept: exact two-factor censuses at GF(13) for witness 2×4² (m ∈ {3,4}), 2×5²
(m ∈ {3,4}), 2×6² (m = 4), 2×3V×witnessB (m = 4), 3×4 MDS pair (m ∈ {5,6});
witness cross-prime probes at $p \in \{7,11,13,17,31\}$ (m = 4); one composite
regrouping instance (m = 4); all subset enumerations complete at the stated
sizes; the four concat anchors. NOT swept: $|S|$ beyond the listed rows;
characteristic-2 fields; infinite fields; factor families beyond those listed;
≥3-factor products beyond the single composite instance; the L-family
counting function; the all-distinct closed form; any asymptotic claim; the
full $d \le 2$ boundary census.
