# Expanded exact series-structure grid

## Fixed ansatz spaces and protocol

[LEMMA] For every displayed ansatz cell with $U$ homogeneous unknown
coefficients, this experiment uses precisely coefficient equations
$0,\ldots,U+1$ as the training prefix and reserves every later safe equation
as a holdout.  Algebraic, Euler, and Mahler matrices use all available stored
orders.  A first-order differential-algebraic row at order $n$ is safe only
through $n=N-2$ when $N$ coefficients of $F$ are known, because the next
coefficient is needed for $[t^n]F'$.

[THEOREM] Let $A$ be one of these finite exact training matrices over
$\mathbb{Q}$ for a fixed displayed ansatz space.  If a stored $U\times U$
minor of $A$ has nonzero determinant, then no nonzero coefficient vector in
that ansatz space annihilates the prescribed training prefix.  Indeed, the
minor gives $\operatorname{rank}_\mathbb{Q} A=U$, so its nullspace is zero.
This theorem is only about the stated finite matrix and ansatz space.

[COMPUTATION] The canonical HT input is `ht_v28.json`, with 29 exact raw
coefficients through $v^{28}$, and the LT input has 33 exact raw coefficients
through $x^{32}$.  The HT prefix through $v^{22}$ agrees exactly with the
frozen predecessor.  The matrix conventions are imported directly from e43
(algebraic and first-order differential-algebraic) and e56 (Euler and Mahler);
no numerical benchmark enters selection or validation.

[COMPUTATION] The evaluated raw-variable grid has
633 admissible cells: HT has
284 and LT has 349.  Its
verdict counts are `{"CANDIDATE_REFUTED_BY_HOLDOUT": 159, "CANDIDATE_SURVIVES": 152, "NO_RELATION_AT_BUDGET": 322}`.
The separately recorded bounded-degree envelopes contain
54519 dimensionally unsupported cells; their
monotone boundary records identify why no training solve was attempted.  Such
skips are not verdicts.

[COMPUTATION] Every `NO_RELATION_AT_BUDGET` row stores a nonzero exact maximal
minor computed after rowwise denominator clearing by fraction-free integer
Bareiss elimination.  Every `CANDIDATE_REFUTED_BY_HOLDOUT` row stores its first
withheld refuting order and exact residual vector on the primitive training
kernel basis.  Exact modular rank was not used to decide any row.  In general
a modular rank can only be a lower bound on rational rank; it would require an
exact confirmation before supporting a no-relation claim.

[COMPUTATION] There are 152 rows labelled
`CANDIDATE_SURVIVES`, but all are certified
`TRUNCATION_UNOBSERVABLE` support kernels, with their complete kernel bases,
residue components, and first testable supported orders in `grid.json`.  Thus
there are 0 unexplained finite-prefix
survivors.  None of the truncation-unobservable rows is called a discovery.

## Re-audit of the e43 LT survivors

[COMPUTATION] The 18 stored e43 LT truncation-unobservable survivors were
replayed in their original $u=x^2$ convention.  The current LT artifact hash
is identical to the one recorded by e43, so no additional equation order has
become available.  All 18
full nullities and zero-column sets agree; the changed-row list is
`[]`.

## Scope

[UNRESOLVED] These are finite exact frontier facts only.  Finite-budget
emptiness is not a non-existence proof for the true infinite free energy, and
a full-rank cell does not prove transcendence, non-algebraicity,
non-differential-algebraicity, or non-D-finiteness.  The calculation neither
claims an exact solution of the simple-cubic 3D Ising model nor excludes
relations outside the displayed ansatz spaces or degree budgets.
