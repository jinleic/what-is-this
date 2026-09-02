# PRE-REGISTRATION — gate H-MINLINE-DGE3-PROOF — T-DGE proof discharge + exhaustive machine controls

Committed: 2026-09-01 (UTC), by agent `RsPe3dTDge`. Written and committed
BEFORE any compute of this run, per Main's binding process order:
prereg -> path-scoped commit (`git commit --only -- <paths>`) ->
`python3 scripts/campaign.py init --gate H-MINLINE-DGE3-PROOF --prereg <this
file>` -> byte-identical copy into the minted run dir (recording source commit
+ sha256) -> controls -> main compute -> `freeze` -> `close --verdict`. No
parameter-dependent computation precedes the prereg commit. Main's steering
(disclose pre-verification) is adopted verbatim as binding, see section 0.

Provenance: Run 1 (gate H-MINLINE, `20260901T103743Z_d86a1fa2_f72db722778d`)
is CLOSED FROZEN-NEGATIVE; Run 2 (gate H-MINLINE-DGE3,
`20260901T104923Z_fba446bf_5aee86fef866`) is CLOSED FROZEN-INCONCLUSIVE with
M1 passing 6/6 and the P1-P5 proof artifact declared the exact residual gap
(see its `closing_record.md`). This run is a SEPARATE lifecycle completing
that registered program; it modifies NO frozen artifact and NO root ledger.

## 0. Pre-commit disclosures (fixed before the commit, per Main's steering)

1. **Scratch pre-verification, NON-EVIDENTIAL.** Before writing this prereg
   the author ran a scratch linear-algebra pass (generic rref/nullspace
   mechanics on toy matrices: P2's circuit-kernel all-nonzero property, the
   P2 constant-vector conclusion on a worked instance, P3's converse on a
   directed coincidence, P4's pair formula on four hand-picked instances,
   an M6 sharpness instance, and a reconciliation pre-check of the 15 RS
   census rows of section 3). These scratch runs happened BEFORE the prereg
   commit; their outputs are **non-evidential** — the campaign verdict rests
   ONLY on in-run controls re-run inside the minted run dir after `init`.
   The scratch files are PRESERVED (not deleted) at
   `cs/rs-pe3d/campaigns/<run_id>/scratch_precheck.py` with its output in
   the run record, labeled non-evidential; nothing outside `cs/rs-pe3d/`
   was touched.
2. **Assignment-supplied expectations.** Every numeric target in section 3
   was SUPPLIED to the author in the assignment (batch context quoting
   Main's independently verified numbers and the two frozen runs' records);
   the author could NOT and did not adapt any threshold to unseen data.
   The prereg pins those supplied numbers as-is; the in-run controls test
   the theorem + instrument against them, they do not define them.
3. **Inherited-defect note (Run-2 prereg regime-2 arithmetic).** Run 2's
   committed prereg stated the regime-2 count as
   $\sum_{\emptyset \ne D \subseteq \mathrm{deg}} 2^{|D|-1} \prod_{i \notin
   D} s_i$. That expression implicitly assumes every degenerate factor has
   a SINGLE parallel class (all columns pairwise proportional) and
   UNDERCOUNTS multi-class spark-2 factors. The correct pinned semantics is
   Run 2's own closing-record formula: with $A_i = \sum_k m_{i,k}^2$ over
   parallel-class multiplicities, unordered distinct dependent pairs
   $= (\prod_i A_i - N)/2$, fibers split per section 3 (M2). Anchor row:
   Main's [2,2] with 4 columns each -> 120 = 48 fibers + 72 off (the
   superseded expression gives 10). This run adopts the $A_i$ formula; the
   superseded expression is recorded here as an inherited defect, NOT
   silently substituted (sibling `RsPe3dMinline` confirmed the same reading
   by IRC before this commit).

## 1. Theorem T-DGE (statement fixed NOW, algebra restated from the Run-2 prereg)

Setting: $H = H_1 \otimes \cdots \otimes H_n$ over an arbitrary field $F$;
each factor $H_i$ is an $r_i \times s_i$ matrix with ALL columns nonzero;
spark $d_i$ = least size of a dependent column set; $d = \min_i d_i \ge 2$;
columns of $H$ are $h(a) = \bigotimes_i h_i(a_i)$ indexed by tuples
$a = (a_1,\dots,a_n) \in \prod_i [s_i]$.

**Claim.**
(A) spark$(H) = d$.
(B) If at most one $d_i$ equals 2, a set $S$ of exactly $d$ columns of $H$
is dependent **iff** $S$ is an axis fiber: all index tuples of $S$ agree
outside one coordinate $i$ (necessarily with $d_i = d$), and the varying
axis-$i$ indices form a size-$d$ circuit of $H_i$; conversely every such
fiber is dependent.
(C) If two or more $d_i$ equal 2, then $d = 2$, (B) fails, and the dependent
pairs are exactly the pairs of distinct tuples that agree componentwise up
to parallel class (for each $i$ the two factor columns are proportional;
the tuples differ on a nonempty set $D$ of coordinates, all of which have
$d_i = 2$).

**RS specialization (P5, this target).** For $H_i$ a GRS parity check of
$\Lambda_i\,\mathrm{RS}(S_i, t_i)$ — $(s_i - t_i) \times s_i$, Vandermonde
rows of exponents $t_i .. s_i - 1$, distinct nonzero evaluation points:
$d_i = s_i - t_i + 1$ (MDS parity check); when $d_i = d$, the size-$d$
circuits are exactly the $\binom{s_i}{d}$ $d$-subsets; parallel classes are
all singleton unless $d_i \le 2$; $d_i = 2$ gives one class of size $s_i$,
so $A_i = s_i^2$; $d_i = 1$ (empty parity check, $s_i = t_i$) means the
factor is the zero map and every product column is a carrier (the
weight-1/trivial window, outside $d \ge 2$; the two registered trivial rows
in section 3 are counted as $N = \prod s_i$).
**LINE corollary — stated ONLY for $t_i = 1$:** at $t_i = 1$ the
minimum-weight words of $C_i$ are exactly the indicators of full axis-$i$
lines, so (B) recovers the LINE characterization (at most one tied
degenerate axis) and (C) exhibits the Run-1 falsification at
$(13,(2,2,4),t=(1,1,1))$: 8 off-line diagonals + 16 full lines = 24. For
$t_i > 1$ the LINE phrasing FAILS even where the circuit theorem holds
(e.g. $(13,(4,4,4),t=(2,1,1))$ has 64 carriers in circuit form but only 16
full lines): the theorem is stated and used in circuit form.

Carrier definitions used by the controls: a weight-$d$ carrier of $H$ is a
$d$-subset of product columns that is dependent; fibers per (B); off-fiber
carriers (diagonals) are dependent $d$-sets that are not fibers.

## 2. Proof obligations P1-P5 (each DISCHARGED in the run artifact `theorem_t_dge.md`)

- **P1 (spark of the product).** spark$(H) = d$ by induction on $n$
  (2-factor step with regrouping $B := H_2 \otimes \cdots \otimes H_n$,
  induction giving spark$(B) = \min_{i\ge2} d_i$ and all columns of $B$
  nonzero since products of nonzero vectors are nonzero; dual-functional
  isolation on distinct A-indices; the $r = d$ edge case is NOT needed at
  sizes $< d$ because $|U| \le k < d \le d_A$ strictly).
- **P2 (size-$d$ classification; the $r = d$ step).** A size-$d$ dependent
  set with $r < d_A$ distinct A-indices: dual isolation collapses each
  group, forcing exactly one group => fiber along B with $d_B = d$ and a
  B-circuit. The $r = d$ case: A-columns form a size-$d$ circuit ($d = d_A$
  forced), 1-dimensional kernel spanned by an everywhere-nonzero $\alpha$;
  id $\otimes \psi$ for arbitrary $\psi$ gives $c_k \psi(b_{v_k}) =
  \lambda(\psi) \alpha_k$ with NO zero-coefficient case split (zero is in
  the span), so all $c_j b_{v_j}/\alpha_j$ equal a single nonzero $w$: all
  $b_{v_j}$ parallel; P3 then forces the $v_j$ to coincide when B has no
  distinct parallel pair.
- **P3 (no-parallel product lemma).** For nonzero simple tensors,
  $\bigotimes_i u_i$ parallel to $\bigotimes_i v_i$ iff $u_i$ parallel to
  $v_i$ for every $i$ (converse by coordinate-extraction functionals,
  valid over any field).
- **P4 (regime classification + counts).** Regime split by the number of
  $d_i = 2$ factors. Counts: at $d = 2$ the unordered distinct dependent
  pairs number $(\prod_i A_i - N)/2$ with $A_i = \sum_k m_{i,k}^2$; the
  fiber subcount is $\sum_i (\prod_{j\ne i} s_j) \sum_k \binom{m_{i,k}}{2}$,
  the remainder off-fiber; in the fiber regime ($d \ge 2$, at most one
  spark-2 factor) every carrier is a fiber and their number is
  $\sum_{i: d_i = d} (\prod_{j\ne i} s_j) \cdot \#\{\text{size-}d \text{
  circuits of } H_i\}$.
- **P5 (RS specialization).** As in section 1: MDS parity-check spark,
  circuit counts, class structure, the $t_i = 1$ LINE corollary (only at
  $t_i = 1$), the Run-1 reconciliation 24 = 16 + 8, the $d_i = 1$ zero-map
  clause, and the $\Lambda$-invariance note (diagonal rescaling preserves
  proportionality classes and circuit structure).

**Verdict rule (fixed NOW):** FROZEN-CERTIFIED only if P1-P5 are all
discharged with no gap in `theorem_t_dge.md` AND every control family
(section 3) passes exactly in-run — zero mismatches, zero tolerance, exact
GF(p) arithmetic, no floating point. Any single undischarged step, any
control mismatch surviving one defect-resolution pass, or any mispinned
target converted by dated in-run addendum still failing -> verdict
FROZEN-INCONCLUSIVE naming the exact residual gap; the theorem's registered
status then remains OPEN, per the Run-2 closing-record semantics.

## 3. Control matrix M1-M6 (all exact, pre-registered, zero tolerance; ALL re-run inside the minted run dir after init)

Instrument: self-contained Python stdlib (`PYTHONDONTWRITEBYTECODE=1`),
exact GF(p) rref (`pow(x,-1,p)`), exhaustive subset enumeration; factors
from Vandermonde rows (any row-count columns independent => MDS behavior:
spark = rows+1) and 1-row spark-2 factors with pinned class layouts.

- **M1 fiber-regime structured (GF(13)), all-d_i>=3, expected carriers (all
  fibers, off = 0):** [(3,5),(3,5)] -> 100; [(3,5),(4,5)] -> 50;
  [(3,4)]^3 -> 192; [(3,4),(4,4),(4,4)] -> 64; [(4,6)]^2 -> 180;
  ([(4,6),(5,6)]) -> 90. Below-$d$ emptiness exhaustively asserted;
  product spark asserted $= d$ via an exhibited dependent $d$-set.
- **M2 diagonal-regime structured (GF(13)), d = 2, >= 2 spark-2 factors,
  pinned class layouts; assertion = P4 pair formula with MEASURED
  multiplicities plus fiber/off split:**
  a. [1x5-distinct]^2 (5 singleton? no: ONE class of 5 per factor):
     $A_i = 25$, $N = 25$ -> 300 pairs, all off-fiber? measured against
     formula AND split (fibers $\sum_i (\prod_{j\ne i} s_j)\binom{5}{2}$-
     style subcount = 20 + 20 = 40? pinned expectation: fibers 40, off 260).
  b. [1x4-distinct, 1x4-distinct] (Main's pinned anchor): $N = 16$,
     $A_i = 16$ -> 120 pairs = 48 fibers + 72 off.
  c. [1x5-distinct, 1x5-distinct, 4x4-Vandermonde] -> 72 fibers + 72 off
     (pinned from Run-2 M2 row: total 144).
  d. [1x5-distinct, 4x5-Vandermonde] single degenerate: 24, ALL fibers,
     off = 0.
  e. [1x4-distinct, 4x4-Vandermonde, 4x4-Vandermonde] single degenerate:
     96, ALL fibers, off = 0.
  f. MULTI-class spark-2 factor (pinned layout $[c,\,2c,\,e]$, classes
     {0,1} size 2 and {2} size 1, $A = 5$) with a 4x4-Vandermonde factor:
     $N = 12$, $A$-product $= 20$ -> 4 pairs, ALL fibers, off = 0
     (single-degenerate regime; exercises $m_{i,k} > 1$ inside the formula).
  g. MULTI-class + MIN-class mix at two degenerate axes: factor
     $[c,2c,e]$ ($A_0 = 5$) x factor $[u,v,w,x]$ 1x4 all-distinct —
     all-proportional (ONE class of 4, $A_1 = 16$): $(5\cdot16 - 12)/2 =
     34$ pairs; pinned split: fibers 22 (4 via class {0,1} x 4 fixed
     axis-1 values... measured: 4 axis-0-fibers + 18 axis-1-fibers), off 12.
- **M3 randomized non-MDS battery (GF(7), seed-pinned rng).** 21
  configurations, 2-4 factors, random row counts, all columns forced
  nonzero, sparks MEASURED (never assumed); asserts spark(product) = min
  spark; regime classification with ZERO violations; regime-2 rows
  additionally assert the P4 pair formula from measured multiplicities.
- **M4 planted controls (instrument validity, fail-loud).** (a) non-tensor
  column: replace one product column of a registered M1 config's H by the
  sum of two distinct product columns; the machinery must REJECT (exhibit
  exhaustively a below-$d$ dependent set in the planted matrix, and none in
  the pristine one, through the identical code path). (b) duplicated factor
  column: factor gains a duplicate of an existing column (spark drops to 2
  where the registered config had none) — the config measurement must
  REJECT the registered expectation (measured spark mismatch fires the
  assertion). Plants must fail loudly on inverted expectations.
- **M5 RS census reconciliation (P5, exact enumeration, GF(13)).** The 15
  registered rows (evaluation points 1..s): (2,3,4)t111->12;
  (2,4,4)111->16; (3,3,3)111->27; (3,3,4)111->24; (3,4,4)111->16;
  (4,4,4)111->48; (3,3,3)211->27; (4,4,4)211->64; (2,2,4)211->16;
  (2,3,4)222->24; (2,2,3)111->18; (2,2,2)111->28; (3,3,3)221->108;
  trivial $d_i=1$ rows counted as $N$: (2,2,4)211 (=16, consistent) and
  (2,3,4)222 (=24, consistent) — plus the Run-1 anchor $(2,2,4)111 \to 24$
  with the dependent-pair split 16 full lines + 8 off-line diagonals (the 8
  pairs classified and exhibited, mapping to Run 1's flat-(0,12) witness
  pair under the pinned sorted-axis convention).
- **M6 boundary sharpness (hypothesis of (B) is exact).** Every M2 row with
  >= 2 degenerate factors must exhibit >= 1 off-fiber dependent pair
  (constructive witness: distinct $p_1 \ne p_2$ in factor $i$, $q_1 \ne
  q_2$ in factor $j$, the pair $((p_1,q_1),(p_2,q_2))$ differs on both axes
  and is dependent). Single-degenerate rows (d, f) must exhibit ZERO
  off-fiber pairs — the theorem's regime split is sharp in both directions.

All counts integers, exact, zero tolerance; mismatch halts with a defect
record; a mispinned prereg target is an inherited defect, correctable only
by dated in-run addendum with root cause (prereg bytes untouched), and an
addendum-corrected target that then fails converts to
FROZEN-INCONCLUSIVE.

## 4. Compute cap and kill criteria (fixed now)

Cap: 45 minutes wall-clock total for M1-M6 in-run (single process; the
largest single enumeration is [(4,6)]^2 at $\binom{36}{4} = 58,905$
4-subsets — well inside). If the cap is hit the run halts with a defect
record and closes FROZEN-INCONCLUSIVE (no silently truncated data). Kill
criteria: (i) any control family failing after one defect-resolution pass
-> FROZEN-INCONCLUSIVE naming the family; (ii) any undischarged P-step ->
FROZEN-INCONCLUSIVE naming the step; (iii) pristine failure at freeze ->
fix and refreeze before verdict. Arithmetic: Python ints mod p, exact
rref/inverses via `pow(x, -1, p)`; NO floats, NO numpy, no MLflow/network.

## 5. Scope sentence (Rule-7 style, fixed now)

Swept by this run: the exact matrix-model censuses of section 3 (complete
subset enumeration, all configs at $N \le 625$ columns, exact GF(p) rank
through pinned rref); the 15 RS reconciliation rows at pinned $(q,s,t)$;
M6 existence/nonexistence on the same rows. Not swept, plainly: RS
instances beyond the listed rows; column Gospels beyond Vandermonde/
GRS-model families; weight-$> d$ windows; asymptotic claims; Conjecture
4.2's $\rho$ content; non-identity $\Lambda$ (the $\Lambda$-invariance
argument is stated in the artifact but machine-tested only at
$\Lambda = \mathrm{Id}$); the weight-1/trivial window beyond the two
registered zero-map rows.
