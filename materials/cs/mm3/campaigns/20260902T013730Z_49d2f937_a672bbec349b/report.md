# Campaign report — gate `mws59-bit-58vs59` (run 20260902T013730Z_49d2f937_a672bbec349b)

## Question (prereg sha `aede812c…`, commit `f1668cd`; amendment 1 commit `b26a5dc`)

Decide the residual bit left by campaign `20260901T133118Z_…7970c5318d5d`: is
`mws59`'s fixed-orientation minimum **58 or 59**?

## Answer — verdict FROZEN-CERTIFIED: the bit is 58

**`mws59`'s fixed orientation has minimum EXACTLY 58 additions (LB = UB = 58),
and an explicit 58-addition circuit is verified end to end.** The published
59-addition scheme is therefore **one addition above the optimum of its own
orientation** — its left stage spends 15 additions where 14 suffice. No
published claim is contradicted: the paper claims 59 additions suffice, and
they do.

    total = C(U) + C(V) + C(output) = 14 + 15 + 29 = 58

| stage | value | how it is established |
|---|---|---|
| C(U) | **14** exact | frozen campaign-1: floor@13 impossible (dual-checker certificate) + a 14-gate circuit extracted from the single admitted auxiliary of 389; re-verified here by exact expansion |
| C(V) | **15** exact | lower: frozen campaign-1 complete scans (aux-1 0/366 at T=13, aux-2 0/79,728 at T=14, dual instrument) ⇒ C(V) >= 15. upper: the paper's printed Table-2 right stage, transcribed verbatim and verified here — 6 u-gates + 9 inline right-operand additions = **15**, reproducing all 23 Table-3 V rows exactly |
| C(output) | **29** exact | lower: C(Wfac) = 15 exact (frozen) + transposition ⇒ >= 29 (preconditions re-audited: 23/23 nonzero rows, rank(Wfac) = 9 over Q). upper: the printed Table-2 output stage — 9 v-gates + 20 combination additions = **29**, reproducing all 9 C-entries exactly |

## The decisive artifact: end-to-end symbolic verification of the assembled circuit

`assembled_58_endtoend.json`. The assembled scheme is

- **left**: the frozen 14-gate U circuit (each gate recomputed from its operand
  structure; all 23 U rows available as ± wires),
- **right**: the printed Table-2 right stage (15 additions),
- **output**: the printed Table-2 output stage (29 additions),

and it was expanded **symbolically over the 81 monomials A_i·B_j**: every one of
the 9 outputs equals its matrix-product bilinear form
`C[i][j] = Σ_l A[3i+l]·B[3l+j]` **exactly** (9/9, zero mismatches), with the
total recounted from the explicit gate lists as 14 + 15 + 29 = **58**. This
check is independent of the DFS/CNF searches, of the factor-block Brent check,
and of any transposition argument.

## Route P — how the two published stages were validated

The paper's printed Table 2 (t0..t6, u0..u5, 23 products with inline signed
operand expressions, v0..v8, 9 output combinations) was transcribed verbatim
from the pinned layout artifact and then:

1. recounted in the fixed model: **left 15 + right 15 + output 29 = 59** —
   matching the paper's claim and the frozen session-3 tally;
2. matched to the Table-3 factor blocks by an exact **product correspondence**:
   coordinate maps (A,B,C) = (id, id, id), a bijection on the 23 products, and
   per-product signs, with **9 products carrying a negated value compensated in
   the output column**; all 621 entry equalities (23·9 + 23·9 + 9·23) hold
   exactly;
3. verified stage by stage (left/right/output all reproduce their Table-3
   targets), then assembled and expanded end to end as above.

## Amendment 1 (disclosed; committed before the amended compute)

The first Route-P attempt **recounted correctly (15/15/29 = 59) but found no
correspondence**, because the pre-registered matcher forced ε·ε′ = +1, forbade a
compensating output-column sign, and assumed both tables share one 3×3
flattening. Those three constraints were unjustified — `(−u)(v)` products are
legitimate with the sign absorbed by the W column, and the frozen record already
documents this convention class (gate A needed the C-permutation
`[0,3,6,1,4,7,2,5,8]`; session 5 recorded a σ² transposition fix). Amendment 1
widened the space to σ_A, σ_B, σ_C ∈ {id, T} plus independent per-product signs,
accepting only on exact equality of every entry, and **strengthened** the R3
wrong-correspondence plant to the widened space. The failed narrow attempt is
preserved verbatim (`verdict_preamendment.json`,
`runner_output_preamendment.json`) and nothing from it is used. In the event the
widening was not even needed for the coordinate maps — the accepted solution
uses (id, id, id); the sole necessary widening was the ε·ε′ = −1 case.

## Route S — registered but not needed

The fallback (aux-3 existential at T = 15) was registered with its closure
argument, an enumeration-only projection of **6.30 × 10^7 ordered chains**
(366 × mean 397.7 × mean 433.2) at a measured **59.2–120 ms/instance**, i.e.
**173–2,100 CPU-hours** — two to three orders of magnitude beyond the 4.0 CPU-h
cap, hence registered as a prefix instrument with a strict
"prefix-statement-only" negative rule plus a positive-only scouting tier.
**Route P settled the bit, so Route S was never entered**: no prefix was run, no
sampling occurred, and no Route-S claim of any kind is made.

## Controls (both directions, all in-run, all passed)

ACCEPT — A1 729/729 Brent over Z in `fmpz` AND pure int; A2 frozen landscape
row reproduced (d = 13/12/14, floors F/F/F, states 132/81/672) and **both
frozen campaign-1 witnesses re-verified by exact expansion** (14-gate U,
15-gate Wfac); A3 extension-positive control on the independent sun56-V target
positive in three instruments (DFS, kissat 4.0.4, CaDiCaL 3.0.1); A4
transposition preconditions re-audited (23/23 nonzero rows per side,
rank(Wfac) = 9 exactly over Q).

REJECT — R1 one-addition-deleted plant on the 14-gate U witness: REJECTED (a
tau class uncovered); R2 operand-sign perturbation: REJECTED by exact
recomputation; **R3 wrong-correspondence plant (two product labels swapped) in
the widened space: REJECTED**; R4 frozen floor@12 on V re-certified INFEASIBLE
with a fresh DRAT→LRAT certificate accepted by both pinned checkers (rc 0/0);
R5 planted wrong d-count assertion fired.

## Independent post-run audit — 6/6 blocks PASS (`postcompute_audit.json`)

Frozen campaign-1 inputs re-read from its verdict artifact; the Table-2 recount
re-derived by an independent tally (15/15/29); the correspondence re-checked as
a bijection with all 621 entry equalities; the end-to-end 58-circuit artifact
re-checked; the certificate replayed through both checkers; verdict consistency
(LB = UB = 58).

## Cost

**2.2 CPU-seconds** of the pre-registered 4.0 CPU-hour cap (Route P is
arithmetic, not search). The expensive part of this result was already paid by
campaign 1's 80,483-instance dual-instrument scans (1.33 CPU-h).

## Frontier statement (bounded, precise)

For the `mws59` factor blocks in their fixed orientation (sigma^0, all-monomial
data triple), in the three-stage linear SLP model (inputs free; gate = x±y;
sign changes and copies free; output counted with the transposition bound):
**the minimum is exactly 58 additions**, attained by the circuit verified here
(frozen 14-gate left stage + the paper's 15-gate right stage + the paper's
29-gate output stage). The published 59 is one above this orientation's
optimum. **58 is not a record** — `paper55`'s 55 stands untouched — and no
claim is made about other mws59 data triples or sigma classes, other
decompositions, non-ternary alphabets, GL(3,ℚ)/GL(3,ℤ) sandwiches, or other
counting conventions.

## Exact fixed-orientation ladder after sessions 9–12

| decomposition | published | certified exact minimum for its own fixed orientation |
|---|---|---|
| `paper55` | 55 | **55** (closed, record) |
| `sun56` | 56 | **56** (closed session 9; published optimal) |
| `mws59` | 59 | **58** (closed here; published is +1) |
| `stapleton60` | 60 | **60** (closed session 11; published optimal) |

## Instrument defects

Amendment 1's cause (my too-narrow registered correspondence) is the only
defect of this run and is disclosed above with its preserved failed attempt.
Carried forward and honoured: no transposition-derived witness is used anywhere
— C(Wfac) = 15 comes from a synthesized circuit and the output-stage 29 from
the printed stage re-expanded and verified.

## exactly-one verdict

**FROZEN-CERTIFIED** — the bit is **58**: `mws59`'s fixed-orientation minimum is
exactly 58 additions (LB = UB = 58), witnessed by an explicit 58-addition
circuit whose 9 outputs were verified equal to the matrix-product bilinear forms
by exact symbolic expansion over all 81 monomials, independent of every search.
