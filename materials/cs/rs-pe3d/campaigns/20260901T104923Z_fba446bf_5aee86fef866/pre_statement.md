# PRE-REGISTRATION — Run 2 — gate H-MINLINE-DGE3 — corrected theorem (proof + exhaustive controls)

Committed: 2026-09-01 (UTC), by agent `RsPe3dMinline`. Written and committed
BEFORE any compute of this run, per Main's binding process order (prereg ->
path-scoped commit -> init -> controls). Run 1 (gate H-MINLINE) is CLOSED
FROZEN-NEGATIVE (campaign `20260901T103743Z_d86a1fa2_f72db722778d`, commit
`8f99430`); this run is a SEPARATE lifecycle with its own verdict, per Main's
verdict-semantics instruction.

## 1. Theorem T-DGE (statement fixed NOW, before controls)

Setting: $H = H_1 \otimes \cdots \otimes H_n$ over any field, each factor
$H_i$ an $r_i \times s_i$ matrix with ALL columns nonzero ("column Gospel":
here the GRS/RS parity or generator matrices), spark
$d_i = \min\{|G| : $ columns indexed by $G$ are dependent$\}$, and
$d = \min_i d_i \ge 2$. A "carrier family" is any $S$ with
$\dim(\mathrm{col}(H) \cap \mathbb F^{S}) > 0$, $|S| = d$ (weight-exactly-$d$
carriers); "fiber" supports are the axis-index sets agreeing outside one
coordinate $i$, with axis-$i$ values a circuit (minimal dependent set) of
$H_i$.

(T-DGE) Exactly one of the following holds:

1. **At most one factor has $d_i = 2$** (at most one nontrivial parallel
   column class in the product): spark$(H) = d$ exactly, and every
   weight-$d$ carrier of $H$ is an axis fiber on a size-$d$ circuit of the
   corresponding factor; conversely every such fiber is a carrier.
2. **Two or more factors have $d_i = 2$**: then $d = 2$ and the weight-2
   carriers are exactly the distinct index pairs agreeing componentwise up
   to parallel class — i.e. $\{p, q\}$ with $p_i = q_i$ off $D$,
   $p_i \ne q_i$ on $D \ne \emptyset$, and for each $i \in D$ the two
   factor columns are proportional; these include genuine off-fiber
   diagonals, and the count is
   $\sum_{\emptyset \ne D \subseteq \mathrm{deg}} 2^{|D|-1}
   \prod_{i \notin D} s_i$ over $\mathrm{deg} = \{i : d_i = 2\}$.

RS specialization (this target): $H_i$ a GRS parity check of
$\Lambda_i \mathrm{RS}(S_i, t_i)$ has spark $d(C_i) = s_i - t_i + 1$ (MDS);
$t_i = 1$ makes optimal words full lines, recovering the LINE form as a
corollary; the falsified H-MINLINE line form fails exactly in regime 2
(tied degenerate axes) — the Run-1 counterexample
$(13,(2,2,4),t=(1,1,1))$: deg $= \{0,1\}$, count
$2^{1}\cdot 4 \cdot$ ... $= 8$ diagonals $+ 16$ lines $= 24$, matching Run 1.

## 2. Proof obligations (human-audited math, machine-anchored)

The proof follows Main's steering route (2-factor lemma + induction), with
audit fixes from Main's independent fable review. Obligations, each to be
stated and discharged in `theorem_t_dge.md`:

- **P1 (2-factor independence below d).** For $A, B$ with sparks
  $d_A, d_B \ge d$ and $k < d = \min(d_A, d_B)$: every $k$
  Kronecker-columns $a_u \otimes b_u$ (all $a_u, b_u \ne 0$) are
  independent. Proof step (per Main): partition the indices by DISTINCT
  A-column; groups $(a, \{u\}, B_u)$ with $|B_u| < d_B$ distinct B columns;
  dual-isolation: for each group pick $\psi$ in theDual$(A)^{\perp}$
  annihilating all other groups' A-columns and separating $a$ (i.e.
  $\psi \in \mathrm{span}$ of rows? audited form: choose functional
  $\psi$ with $\psi(a) = 1$, $\psi(\text{other groups' A columns}) = 0$,
  possible since the A columns of other groups plus $a$ are independent
  when their total is $< d_A$); apply id $\otimes \psi$ to the relation:
  $\sum_{u \in B_u\text{-group}} c_u \psi(a_u) b_u = 0$ collapses the whole
  relation to one group, where the B columns are distinct and $< d_B$
  many, hence independent, forcing $c_u = 0$. Induction over groups.
  **Audit point (Main's fable):** the r = d edge case: when a size-d
  dependency has $r = d_A$ distinct A indices, the isolation functional may
  not exist; handled by the circuit lemma P2 instead. Zero coefficients:
  $\psi$ may annihilate the chosen $a$ too (then the group contributes
  $\lambda(\psi) = 0$ and is merely dropped) — the separation argument only
  needs SOME group surviving with all $c_u b_u$ visible.
- **P2 (size-d dependency structure, r = d case).** A size-$d$ dependent
  set of Kronecker columns with $r = d = d_A$ distinct A columns: the
  A-part is itself a size-$d$ circuit of $A$ with 1-D kernel
  $\mathrm{span}(\alpha)$, $\alpha$ everywhere nonzero; apply
  id $\otimes \psi$ for arbitrary $\psi$: $(c_k \psi(b_k))_k \in
  \mathrm{span}(\alpha)$, so $c_k b_k / [\alpha_k \ne 0]$ is constant in
  $k$ — all B columns proportional. Spark $\ge 3$ of $B$ (or: $d_B > 2$)
  forbids DISTINCT proportional columns, so the B index is unique, and the
  dependency is an A-circuit on a fixed B index: a fiber. Ditto with A/B
  roles exchanged. **Audit point:** this yields the fiber classification
  whenever BOTH factors have $d_i \ge 3$ OR exactly one has $d_i = 2$
  (that one being the "B-side" whose parallelism is absorbed: with a
  single spark-2 factor the size-2 dependency is the fiber on that axis
  and no diagonal arises — verified: $(13,(2,3,4))$ 12/12,
  Main's $[2,4] \to 24/0$, $[2,4,4] \to 96/0$).
- **P3 (no parallel columns from non-parallel factors).** If each $H_i$
  has no two distinct proportional columns, then so does any product
  $H \otimes H'$ — a Kronecker-column pair proportional forces each
  factor pair proportional. This is what makes $d_i \ge 3 \,\forall i$
  (no factor fails) the clean regime.
- **P4 (regime 2 classification + count).** With $\ge 2$ factors having
  spark 2: $d = 2$; a pair of product columns is dependent iff each factor
  pair is proportional (P1 contrapositive at $k = 2$); agreeing up to
  parallel class per coordinate and distinct overall. Count: for each
  nonempty $D \subseteq \mathrm{deg}$ of axes where the pair DIFFERS
  (necessarily along parallel classes), $2^{|D| - 1}$ orientation
  identifications times $\prod_{i \notin D} s_i$ fixed coordinates.
- **P5 (RS specialization).** GRS parity-check spark $= s_i - t_i + 1$;
  $t_i = 1$ line corollary; the falsified H-MINLINE failure localized to
  regime 2 (Run-1 instance reconciled: 24 = 16 + 8).

Verdict rule (fixed now): T-DGE is promoted only if P1-P5 are all discharged
with no gap in the written artifact; the exhaustive matrix-model controls
below must pass exactly. Any gap: theorem stays OPEN with the missing lemma
stated precisely, and the campaign closes FROZEN-INCONCLUSIVE.

## 3. Exhaustive controls (exact GF(p), fixed now)

Model instrument (Main's protocol, re-derived): for random/exact matrices
over GF(p), measure spark per factor (exhaustive), spark of product,
enumerate ALL weight-$d$ subsets of product columns (exact rank), classify
fiber vs off-fiber, and compare against T-DGE's closed forms.

- **M1 fiber-regime structured:** $d_i \ge 3$ configs
  $[3,3]$ (N=25, expect dep 100 off 0), $[3,4]$ (50, 0), $[3,3,3]$ (192, 0),
  $[3,4,4]$ (64, 0), $[4,4]$ (180, 0), $[4,5]$ (90, 0) — Vandermonde
  factors over GF(13), ALL subsets of size $d$ enumerated.
- **M2 diagonal-regime structured (d=2):** $[2,2]$ (expect 120 = 48 fibers
  + 72 off), $[2,2,4]$ (144 = 72 + 72), $[2,4]$ (24, 0),
  $[2,4,4]$ (96, 0).
- **M3 randomized non-MDS:** 21 configs over GF(7), 2-4 factors, mixed
  measured sparks (including spark-2 with a single degenerate axis and
  none): check spark(product) = min spark exactly and T-DGE regime
  classification with ZERO violations.
- **M4 RS-P3d census reconciliation:** the Run-1 instances
  $(13,(2,2,4))$, $(13,(2,3,4))$, $(13,(3,3,4))$ at $t = (1,1,1)$ through
  the T-DGE count formula (fiber and diagonal regimes): 24, 12, 24
  respectively; plus Main's additional RS rows
  $(2,4,4) \to 16$, $(4,4,4) \to 48$, $(3,3,3) \to 27$,
  $(3,3,4) \to 24$ (all t=1), and $t \ne (1,1,1)$ circuit-form rows
  $(3,3,3),t=(2,1,1) \to 27$, $(4,4,4),t=(2,1,1) \to 64$,
  $(4,4,4),t=(2,1,1)$ line-form comparison ($\ne$ 16 full lines — line
  phrasing fails when $t_i > 1$ even where the circuit theorem holds),
  $(2,2,4),t=(2,1,1) \to 16$, $(2,3,4),t=(2,2,2) \to 24$,
  $(2,2,3) \to 18$ (t=1), $(2,2,2) \to 28$ (t=1), $(3,3,3),t=(2,2,1) \to 108$.
- **M5 planted objects (instrument controls):** pure-tensor column ACCEPT;
  non-tensor sum column REJECT; duplicated factor column REJECT (plant
  off-class: produces below-$d$ dependency, exhaustively visible);
  fail-loud on any expectation mismatch.
- **M6 prior-instance cross-check:** Run-1's exact numbers (24/16/8 with
  the 8 diagonals listed) must re-derive from the T-DGE classification of
  the RS instance — same objects, new description.

All counts integers, zero tolerance; mismatch = halting defect record.

## 4. Rule-7 scope sentence (fixed now)

Swept: the matrix-model censuses of section 3 exhaustively (complete
subset enumeration per instance at $N \le 625$ columns... instances sum to
< 2000 subset rank tests per config, all exact GF(p) rank by rref); the RS
count reconciliations listed (closed forms + exact enumeration where
$N$ small). Not swept, plainly: numeric-range census of every RS instance
beyond the listed rows; non-GRS column Gospels beyond the model families;
asymptotic claims; Conjecture 4.2's $\rho$ content; weight-$> d$ windows;
non-identity $\Lambda$ beyond rank argument (P5 notes $\Lambda$ preserves
proportionality classes: diagonal rescaling maps columns to columns and
proportional pairs to proportional pairs — stated in the artifact).
