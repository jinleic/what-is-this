# Gate-B encoding audit (ordered by Main, 2026-08-30) — verdict: SOUND, both directions

Trigger: arXiv:2607.29291 (Palladinos) found that ten Heule–Kauers–Seidl
"expected-UNSAT" Challenge-2 formulas are in fact SATISFIABLE, because
*"their hardcoded type-3 pairings are imposed by positive unit clauses on the
621 base variables: the formulas require selected incidences but do not forbid
additional type-3 incidences."* Main's instruction: audit gate B's own encoding
for the same defect shape — REQUIRE without FORBID — without re-running gate B.

Audited object: `cs/mm3/src/gate_b_encodings.py`, functions `build_floor_cnf`
(the CNF that produced the CaDiCaL/kissat UNSATs certified by drat-trim +
lrat-check in campaign `2026-08-30T031544Z_e0f3f117_c9df97a8bf3a`) and
`build_floor_ilp_data` / `solve_floor_ilp` (the HiGHS route).

Variables: `Y[g,c]` = "slot g produces needed class c" (g < T, c < n = d(F));
`P[g,c,i]` = "slot g produces c via representation i", where representation i
ranges over `reps[c]` = **every** unordered pair of base classes (a,b) with
±a ± b = ±c, enumerated exhaustively by `prep()`.

## (a) Which clauses FORBID and which merely REQUIRE

| family | clause | role |
|---|---|---|
| (1) | `¬Y[g,c₁] ∨ ¬Y[g,c₂]` for all c₁<c₂, all g | **FORBIDS** a slot producing two classes |
| (2) | `⋁_g Y[g,c]` for each needed c | REQUIRES coverage of every needed class |
| (3a) | `¬Y[g,c] ∨ ⋁_i P[g,c,i]` | **FORBIDS** production without a representation |
| (3b) | `¬P[g,c,i] ∨ Y[g,c]` | **FORBIDS** a representation flag without its production (full reification, both directions present) |
| (4) | `¬P[g,c,i] ∨ ⋁_{h<g} Y[h, cls(operand)]` for each non-input operand | **FORBIDS** using a representation whose operand class was not produced in an *earlier* slot — this is the "…and only these" family whose absence broke Challenge-2 |
| units | `¬Y[g,c]` when `reps[c] = ∅` (CNF); `y ≤ 0`, `p ≤ 0` rows (ILP) | **FORBID** (negative units only) |

**The Challenge-2 defect shape is absent by construction.** That defect was
*positive* unit clauses asserting selected incidences with no companion clause
excluding extra ones. This encoding contains **no positive unit clauses at all**;
its only units are negative (forbidding unreachable classes). The single REQUIRE
family, (2), needs no companion FORBID: slots are capped at one class each by
(1) and there are exactly T = d(F) slots for d(F) classes, so "more than
required" is not expressible.

## (b) Two-direction correspondence (both stated, both hold)

*Forward — model ⇒ circuit.* Given a satisfying assignment: by (1) each slot g
carries at most one class c; by (3a) some `P[g,c,i]` is true, giving a
representation (a,b) with ±a ± b = ±c; by (4) every non-input operand class was
produced at some slot h < g. Emit gates in slot order, gate g := ±a ± b. Inputs
are free and a produced class makes both signs available at zero cost (sign
changes are free in the model), so each gate's operands are available when
needed. By (2) every needed class is produced, and targets that are ± inputs
need no gate. Hence a circuit with ≤ T gates computing all 23 targets exists.

*Reverse — circuit ⇒ model.* This is the direction the Challenge-2 authors left
unstated. Let a T-gate circuit exist with T = d(F). By the floor lemma each gate
creates at most one new needed direction class, and d(F) classes must be
created, so with exactly d(F) gates the map "gate ↦ class it creates" is a
bijection and **every gate value is ± a needed target class** — not an
assumption but a consequence. Topologically order the gates as slots, set
`Y[g,c] = 1` for the class created at slot g (satisfies (1) and (2)). Each gate
is ±a ± b of available values; the unordered pair (canon a, canon b) is in
`reps[c]` because `prep()` enumerates **all** such pairs — so the matching
`P[g,c,i]` can be set true, satisfying (3); and its non-input operands were
produced strictly earlier, satisfying (4). So every floor circuit has a model.

Both directions therefore hold, and UNSAT ⇔ no d(F)-gate circuit exists. The
certified UNSATs decide exactly the intended question.

**Scope note (unchanged, and load-bearing).** The encoding is used **only at
T = d(F)**, where the value restriction is a theorem. At T = d+1 the same
encoding is deliberately *too strong* (it forbids auxiliary values), which is
precisely what campaign `…031544Z…`'s control observed: raw floor CNF at d+1 is
UNSAT, while the aux-1 extension at d+1 is SAT and admits the paper's own
witnesses. That control is the encoding-adequacy evidence, and it is consistent
with this audit rather than a substitute for it.

## (c) Under-constraint found? **No.** Nothing to escalate.

The one modelling premise shared by all three procedures (DFS, ILP, SAT) is
"availability is at class level because sign changes are free". That premise is
the paper's own cost model (negation free, established in the target's
convention), not an artifact of any encoding — and Main's warning stands
independently: three solvers agreeing is one premise checked three times. The
premise here is justified by definition of the model, and the reverse-direction
argument above is what connects the model to actual circuits.

## (d) Solver-tolerance exposure (Main's second question, 2026-08-30)

- **Gate C (this campaign) uses no solver and no floating point anywhere.**
  All arithmetic is exact: Python arbitrary-precision integers in the sweep and
  the DFS, `python-flint` `fmpz` in `verify_anchors.py`. Factor entries are in
  {−1,0,1} and Brent sums are bounded by 23 in absolute value. There is no
  tolerance to misconfigure and no float to lose.
- **Gate B's HiGHS route:** every variable bound is exactly `[0.0, 1.0]`, every
  matrix coefficient is `±1`, every RHS is `0` or `1` (see
  `build_floor_ilp_data` / `solve_floor_ilp`). The smallest nonzero magnitude in
  the model is 1.0, seven orders of magnitude above HiGHS's default 1e-7 primal
  feasibility tolerance, so the failure mode Main reported (tolerance equal to
  the quantity being constrained) cannot arise here.
- Moreover the HiGHS verdict is **not load-bearing**: `solve_floor_ilp` skips
  the rows whose `rhs is None`, which *relaxes* the model, and a relaxation
  being INFEASIBLE still implies infeasibility of the tighter model — the sound
  direction. The same three floor questions additionally carry binary DRAT
  proofs from kissat 4.0.4 replayed by two independent checkers (`drat-trim`,
  `lrat-check`) plus the complete memoized DFS. A tolerance artifact in HiGHS
  could not flip gate B's conclusion.

**Verdict: gate B is strengthened, not merely re-asserted.** C(U)=13, C(V)=14,
C(W-factor)=14 stand, and the encoding is now audited in both directions with
the Challenge-2 defect shape explicitly excluded.
