# Amendment 1 — control R5's instance was mis-specified; replaced by a genuinely known-SAT instance

**Created 2026-09-02, AFTER the first (aborted) run and BEFORE any amended
compute.** Committed path-scoped before the amended runner executes.
Campaign: `20260902T023856Z_28ff8823_43330ec64df9`, gate
`laderman23-ladder-total`, prereg sha
`2af3b859eff465c3146cda7a9b7cc4ef883137267cda7eced132449fdf4ca4c2`.

## What happened

The first run passed every anchor and the accept-side controls and then aborted
at **control R5**, before the censuses and before any claim. **No verdict was
produced.** The aborted output is preserved verbatim as
`verdict_preamendment_aborted_R5.json` and the pre-amendment runner as
`laderman_ladder_run.py.preamendment`. Everything the aborted run did establish
is unchanged and is re-derived in the amended run:

- triple sha `522ba07f4f1784ac…`, Brent **729/729** int and fmpz;
- `d(U) = d(V) = d(WFac) = 14`, all three floors impossible, 2401 DFS states
  each;
- preconditions: 23/23 nonzero rows per side, `rank(WFac) = 9`;
- **Laderman's printed basic form recounts to `28 + 28 + 42 = 98` additions**;
- A2 the known-true `paper55` printed circuits accepted at their exact counts
  (13 and 14 gates, recount `13 + 14 + 28 = 55`);
- A3 frozen rows reproduced `(12,False) (13,False) (13,False)`;
- A4 relabeling invariance holds.

## The defect: my R5 instance is not a known-SAT instance in this encoding

Pre-statement §11 specified R5 as *"an instance known SAT (a side's classes at
T = 15 + 5 free extra slots, trivially schedulable) must come back SAT from
BOTH instruments"*. That premise is **false for this encoding**, and the
encoder was right to say UNSAT.

Reason: the slot-availability CNF's operand pool is `base = inputs ∪ tau`, i.e.
**only free inputs and needed classes may be used as gate operands** — that is
exactly the closure property that makes the floor question decidable. Adding
slack slots therefore adds no expressive power: whether all of `tau` can be
created from `inputs ∪ tau` is independent of `T` once `T >= |tau|`. So `U`'s
14 classes at `T = 19` are UNSAT precisely because they are UNSAT at `T = 14`
(the floor verdict the DFS already reports). My "trivially schedulable" claim
confused *slots* with *available operands*.

This is a **control-design defect in my pre-registration**, caught by the
control itself before any claim — not a defect of the encoder, and not a fact
about Laderman's scheme.

## The replacement (fixed before the amended run)

R5's purpose is unchanged: guard against an encoder that answers UNSAT
unconditionally. The instance is replaced by one that is **known SAT by
construction**:

> **R5 (amended).** Build a *creatable chain* of 7 classes over the 9 free
> inputs: `c_1 = canon(e_0 + e_1)`, and `c_{i+1} = canon(c_i + e_{i+1})` for
> `i = 1..6`. By construction each `c_i` is a signed sum of two quantities
> available before it, so `tau = {c_1..c_7}` is schedulable in exactly 7 slots.
> **Both instruments must return feasible** at `T = 7`: the CNF must be SAT and
> the complete subset-DFS must report `floor_achievable = True`. Either
> instrument answering infeasible aborts the campaign.

This is strictly stronger than the original intent, because it exercises **both**
instruments in the positive direction on the same instance (the original only
demanded SAT from the CNF).

**Calibration disclosure.** Before writing this amendment I ran exactly this
chain instance through both instruments to confirm the replacement control is
well-posed: CNF **SAT**, DFS **achievable = True**, 359 vars / 695 clauses.
That is control calibration on a synthetic object; it touches none of the
target's censuses, universes, or counts, and it is re-run in-run as the actual
control.

Everything else in the pre-statement stands unchanged: the object, the model
and addition-count semantics, the stagewise reduction, the closure argument,
the hash-pinned universes (3 x 408), both instruments and their versions, the
mandatory per-instance agreement, the cross-solver control, the certificate
policy and pinned checkers, the symmetry reductions, the synthesis seeds
(master `20260902`, 240 restarts/side), the 2.0 CPU-h budget, the checkpointing
scheme, controls A1-A6 / R1-R4, the verdict rule, and §13's list of what may
never be claimed.
