# PREREG-RSPE3D-6-H2-CLOSURE — gate **H-MIX-H2-CLOSURE**

Committed: 2026-09-01 (UTC), agent `RsPe3dH2`. This file is written and
path-scoped committed before any parameter-dependent computation for this run.
The binding order is: preregistration → `git commit --only` of this file →
`campaign.py init --gate H-MIX-H2-CLOSURE` → byte-identical preregistration copy
in the minted run directory with source commit and SHA-256 → controls and
compute → proof record → `freeze` → `close --verdict` with exactly one verdict.
No frozen campaign directory is modified.

Lineage: the frozen run
`20260901T132343Z_4f169cbb_2f81dd8f2d8d` (gate `H-MIX-CROSSING`) proved Theorem
X's forward crossing construction and its converse except one general residual,
Hypothesis H2. The target-level `SCOPE_NOTE_H_MIX_CROSSING.md` and `state.json`
register that residual. The five-step argument supplied to this agent by owner
steering is a candidate, not authority: every implication below must be checked
independently in `theorem_h2_closure.md`. If any implication fails, the run names
the first exact gap and H2 remains open.

## 1. Setting and fixed H2 statement

Let $A$ and $B$ be matrices over a field, with all columns nonzero and
$d_A=d_B=3$. Let $H=A\otimes B$, $m_*=d_A+d_B-2=4$. Use the frozen kernel
identity

$$\ker(A\otimes B)=(\ker A\otimes F^{s_B})+(F^{s_A}\otimes\ker B).$$

Thus a relation matrix can be written $\Gamma=C_1+C_2$, where each nonzero
column of $C_1$ lies in $\ker A$ and each nonzero row of $C_2$ lies in
$\ker B$. For a size-four circuit $S=\operatorname{supp}\Gamma$ of profile
$(3,3)$, let $k_1$ be the number of nonzero columns of $C_1$ and $k_2$ the
number of nonzero rows of $C_2$.

**Hypothesis H2 (fixed statement).** If such a circuit relation admits a
decomposition with $(k_1,k_2)=(2,2)$, then $S$ is a crossing support

$$S(T,w,Z,x)=\big((T\setminus\{w\})\times\{x\}\big)\cup
\big(\{w\}\times(Z\setminus\{x\})\big),$$

where $T$ is a three-circuit of $A$, $Z$ is a three-circuit of $B$, $w\in T$,
and $x\in Z$. Equivalently, the $(2,2)$ decomposition corner produces no new
size-four profile-$(3,3)$ circuit support, and its relation is proportional to
the crossing relation.

This run proves or leaves open exactly H2. It does not re-prove frozen Lemmas P,
Q, Q2, Q3, or R except where their statements are needed to state the
consequence.

## 2. Five analytic obligations

Write the two active $C_2$ rows as $U=\{u_1,u_2\}$ and the two active $C_1$
columns as $V=\{v_1,v_2\}$, with distinct indices in each set. Put
$X=U\times V$ (four cells), $P=\operatorname{supp}C_1$,
$Q=\operatorname{supp}C_2$, $O=P\cap Q$, $D_1=P\setminus Q$, and
$D_2=Q\setminus P$.

1. **Confinement, including noncancelling overlap inside $S$.** Prove separately
   that every overlap cell, whether inside or outside $S$, lies in $X$ because
   its row and column are active; and that every supported cell outside $S$
   belongs to both summands and cancels there. This must justify $O\subseteq X$
   without assuming that an inside-$S$ overlap cancels.
2. **Tight count.** Each active $C_1$ column has support at least three and each
   active $C_2$ row has support at least three, so $|P|,|Q|\ge6$. Prove that
   $D_1$ and $D_2$ are disjoint by definition and are genuinely subsets of
   $S$, then use $|O|\le4$ to obtain
   $$4=|S|\ge |D_1|+|D_2|=|P|+|Q|-2|O|\ge4.$$
3. **Equality structure.** Analyze equality in every link of that chain. It
   must force $|P|=|Q|=6$, $|O|=4$, $O=X$, and
   $S=D_1\mathbin{\dot\cup}D_2$. In particular $X$ must lie in both supports,
   not merely bound their intersection, and $X\cap S=\varnothing$, so all four
   cells of $X$ actually cancel. Deduce that the two column supports are
   $T_i=\{u_1,u_2,w_i\}$ and the two row supports are
   $Z_i=\{v_1,v_2,x_i\}$, each of size three.
4. **Profile collapse, with exclusions explicit.** Prove
   $w_i\notin\{u_1,u_2\}$ and $x_i\notin\{v_1,v_2\}$ from the size-three
   support statements. Then prove profile $(3,3)$ forces exactly
   $w_1=w_2=:w$ and $x_1=x_2=:x$. Explicitly rule out purported alternatives
   such as $w_1=u_2$ and their symmetric column versions. Deduce
   $$S=\{(w,v_1),(w,v_2),(u_1,x),(u_2,x)\}.$$
5. **Circuit and crossing conclusion.** Prove that a support-three nonzero
   vector in $\ker A$ is a circuit support, not merely a dependent set: it is
   dependent and every proper subset is independent by $\operatorname{spark}
   (A)=3$. Do the same for $B$. Therefore
   $T=\{u_1,u_2,w\}\in\mathrm{Circ}_A(3)$ and
   $Z=\{v_1,v_2,x\}\in\mathrm{Circ}_B(3)$ and the displayed $S$ is precisely
   $S(T,w,Z,x)$. Finally use the one-dimensional relation space of a circuit to
   prove proportionality to the crossing relation.

**Promotion rule.** The verdict may be `FROZEN-CERTIFIED` only if all five steps
are discharged in writing with no gap, including every emphasized subpoint.
If any step fails, the only admissible substantive verdict is
`FROZEN-INCONCLUSIVE`, naming that step and leaving H2 open regardless of
numerical evidence.

## 3. Exact in-run control matrix

All controls run after `init`, in exact integer arithmetic modulo the stated
prime, standard library only. Every positive path has a known-true ACCEPT and
the run contains planted REJECTs.

- **T1 — owner GF(13) ACCEPT.** Use $A=B=[(1,t)]_{t=1}^4$ as $2\times4$
  Vandermonde factors; $T=Z=\{0,1,2\}$; $u_1=0,u_2=1,w=2$;
  $v_1=0,v_2=1,x=2$; and $c=\delta=(1,11,1)$. Set
  $\alpha=(\delta_{v_1},\delta_{v_2})$ on columns $v_1,v_2$ and
  $\beta=(-c_{u_1},-c_{u_2})$ on rows $u_1,u_2$. Assert the exact support is
  `[(0,2),(1,2),(2,0),(2,1)]`, equals $S(T,w,Z,x)$ as a set, satisfies
  $A\Gamma B^{\mathsf T}=0$ and the coefficient-sum product relation exactly,
  and is a genuine circuit (rank three and every three-subset independent).
- **T2 — exhaustive $(2,2)$ parameter sweep.** At each
  $p\in\{7,11,13\}$ sweep two distinct factor pairs:
  `V4xV4`, with both factors $[(1,t)]_{t=1}^4$, and `V3xV5`, with point sets
  $(0,1,3)$ and $(0,1,2,4,5)$. For every unordered active row pair $U$, every
  unordered active column pair $V$, every ordered choice of third support
  indices $(w_1,w_2)\in([s_A]\setminus U)^2$, every ordered
  $(x_1,x_2)\in([s_B]\setminus V)^2$, use the unique normalized circuit
  direction on each $U\cup\{w_j\}$ and $V\cup\{x_i\}$. Sweep every nonzero
  scalar quadruple modulo its common global scalar by fixing the first column
  scalar to one and enumerating the other three over $F_p^\times$. This is a
  bijective normalization because every component is nonzero; it represents
  every raw scalar quadruple exactly once modulo global scaling. The exact
  normalized tuple totals are pinned by the enumeration formula:
  `V4xV4`: $576(p-1)^3$; `V3xV5`: $270(p-1)^3$.

  Retain every tuple for which all four cells of $X$ cancel. For every retained
  support of size four and profile $(3,3)$, assert genuine circuitness and
  **set membership**, not count-only agreement, in the independently
  constructed crossing family. Assert the set difference “retained
  profile-$(3,3)$ supports minus crossing supports” is empty and that the
  distinct retained population equals the full crossing family (144 for
  `V4xV4`, 90 for `V3xV5`). Record all other retained profiles without
  promoting them. This is exhaustive for H2 because analytic Step 3 proves
  that every H2 candidate has exactly these four size-three component supports
  and four nonzero component scalars.
- **T3 — previous-campaign census rerun.** At GF(13), reuse by value (no import
  from the frozen directory) the previous witness
  $A=((1,1,1,1),(1,2,3,4))$ and
  $B=((1,-2,1,0),(0,1,-2,1))$. Exhaust all $\binom{16}{4}=1820$ four-subsets,
  identify circuits by exact rank and all proper-subset checks, and assert the
  profile-$(3,3)$ circuit set equals the independently constructed crossing
  family: 144 on each side and zero measured-minus-constructed residual. Pin
  the prior total 156 and profile-$(4,4)$ residual 12.
- **T4 — planted REJECTs.** (a) Perturb one nonzero coefficient of T1's forced
  $\Gamma$ by $+1$ and assert the perturbed coefficient vector no longer gives
  $A\Gamma B^{\mathsf T}=0$ (the support itself remains the known circuit).
  (b) Duplicate one factor column into a fifth column and assert exact spark
  measurement drops from three to two and a dependent pair is found.
- **T5 — fail-loud tensor guard.** Independently form expected row-major
  Kronecker coordinates and compare them with the instrument's product columns
  for factors of ambient row dimensions two and three. Assert product ambient
  dimension $2\cdot3=6$, while the planted block-concatenation path has
  dimension $2+3=5$ and different coordinates. This control fails immediately
  if block concatenation is silently substituted for Kronecker product.

## 4. Arithmetic, budgets, defects, and kill criteria

- Exact arithmetic only: Python integers reduced modulo $p$, Gauss–Jordan rank
  with `pow(pivot,-1,p)`, no floats in mathematical decisions, no NumPy or
  external algebra package. `sys.dont_write_bytecode=True` is set before any
  possible frozen-directory import; the run imports nothing from a frozen
  directory. Execution sets `PYTHONDONTWRITEBYTECODE=1`.
- One process at `nice -n 15`, with
  `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=VECLIB_MAXIMUM_THREADS=1`.
- Hard budget: **600 CPU seconds** and **900 wall-clock seconds**. The instrument
  checks both during every exhaustive family and fails before promotion if
  either cap is reached. The scalar-sweep totals above are fixed; no adaptive
  parameter extension is allowed.
- Any defect is disclosed in-run; broken output is preserved; every affected
  family is rerun from its beginning after the repair. No failed assertion is
  weakened or silently repinned.
- Kill criteria: any analytic gap; any ACCEPT or REJECT failure; any nonempty
  set-level residual; failure to enumerate the pinned tuple total; arithmetic,
  tensor-layout, bytecode, thread, CPU, or wall-budget violation; or unresolved
  implementation defect. On kill, preserve the strongest completed evidence
  and close `FROZEN-INCONCLUSIVE` naming the first load-bearing failure.

## 5. Consequence and nonclaims fixed in advance

If H2 is discharged, combining it with frozen Lemmas P, Q, Q2, and R completes
Theorem X's converse for every $d_A,d_B\ge3$ under Theorem X's standing
hypotheses: the $(3,3)$, $(2,2)$-decomposition corner adds no noncrossing
support. It does **not** settle the spark-two boundary branches, an all-distinct
closed form, or general results for three or more factors. Those remain open
and must be stated explicitly in the theorem record and, after a certified
close only, in the target-level scope note and state record.
