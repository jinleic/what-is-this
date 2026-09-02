# Campaign report — gate `sun56-gap-total55` (run 20260901T123618Z_3f41e268_0264d9d39939)

## Question (pre-registered: pre_statement.md, sha 7854d469…)

Does a 55-addition three-stage circuit exist for the `sun56` decomposition in
Sun's fixed orientation? Frozen prior state: sun56 = LB 55 / UB 56, with the
all-monomial data triple the unique LB-55 data triple of the frozen sun56 census.

## Answer

**No.** The 55-move does not exist for this orientation. Reduced exactly:

    total = C(U) + C(V) + C(Wfac) + 14 = 13 + C(V) + 30

with C(U) = 13 (frozen checker-certified floor@12 impossibility + Sun's
13-gate witness, exactly re-expanded), C(Wfac) = 16 (frozen achievable floor,
replayed) and output >= 30 (transposition, precondition audited: 23/23 active
products, exact rank 9). The whole gap is C(V) ∈ {12, 13}.

C(V) = 13:

- floor@11 impossible (frozen DRAT/LRAT certificates, both checkers; re-verified
  in-run R3 with fresh certificates DRAT→LRAT, drat-trim `s VERIFIED` + lrat-check rc 0);
- C(V) <= 12 ⟺ some single auxiliary direction a ∉ tau(V) ∪ inputs admits a
  12-slot schedule over tau(V) ∪ {a} (closure argument in pre-statement §4:
  12 gates covering 11 tau classes have exactly one non-tau gate value, whose
  class must be a signed sum of two elements of inputs ∪ tau(V)); the
  auxiliary universe is exactly 338 elements (CSV sha256 c10bdeea…), computed
  by two independent enumeration loops;
- all 338 instances decided INFEASIBLE by BOTH instruments, agreeing on every
  instance: (1) the complete memoized subset-DFS (gate_b_floor.subset_dfs) at
  T = 12; (2) a freshly written slot-availability CNF solved by kissat 4.0.4.
  Checkpointed per-instance in scan_checkpoint.jsonl (338/338 rows);
- Sun's own 13-gate V-circuit is the witness: expanded exactly, reproduces all
  23 V rows, 13 distinct gate classes = 11 tau + 2 auxiliary — consistent with
  the two-aux floor-achievable control (T = 13 SAT in both instruments).

Hence sun56 fixed-orientation minimum = 13 + 13 + 16 + 14 = **56 exactly**:
Sun's published 56-addition scheme is optimal for its own orientation, and the
gap sun56 = LB 55 / UB 56 closes upward to LB = UB = 56.

## Frontier statement (bounded, precise)

`sun56` (sigma^0 monomial orientation of Sun's factor blocks, three-stage
linear SLP model, additions incl. subtractions, inputs free, sign changes
free, output via transposition accounting): **no circuit with <= 55 additions
exists**. This is NOT a universal claim: other orientations of the same
tensor, other decompositions, other models are untouched. Combined with the
frozen census (campaign 192400Z: all non-monomial sun56 data triples >= 56),
this eliminates the last candidate for a 55-total inside that frozen swept
set — the corollary "the swept set contains no 55-addition scheme other than
the published paper55 one" now holds with sun56 closed at 56.

## Search space, budget, measured cost (assignment points 2-3)

- Space: the exactly enumerated 338-element auxiliary universe (hash-pinned
  in the pre-statement BEFORE compute) x {12-slot floor-lease question},
  plus the already-frozen floor@11 and U/Wfac side decisions.
- Dual instrument per instance with mandatory agreement (disagreement ⇒ abort).
- Budget: pre-registered 3.0 CPU-h hard cap; measured 73.5 s total
  (well under), single-thread, nice -n 15, PYTHONDONTWRITEBYTECODE=1.
- Checkpointing: scan_checkpoint.jsonl, one fsynced JSON line per instance;
  the run is resume-safe and the evidence is the complete 338/338 prefix.

## Controls (assignment point 4) — all in-run, both directions

ACCEPT: A1 anchors (729/729 Brent over Z in fmpz AND pure-int; split recount
13/13/30 = 56); A2 witnesses (U and V circuits reproduce all rows by exact
expansion; counts recomputed 13/13/30); A3 extension-positive (tau ∪ {Sun's
two aux} schedulable at T = 13 in both instruments — the witness-shape
instance); A4 frozen-row reproduction (d/floor/states 12-33, 11-9, 16-17 and
total_lb 55 re-derived from raw factors through the same code path).

REJECT: R1 one-addition-deleted plant (V gate 6 deleted ⇒ verifier REJECT by
exact expansion); R2 operand-perturbed plant (V gate 3 operand-1 sign flipped
⇒ verifier REJECT); R3 floor@11 must be infeasible in both instruments with
fresh dual-checker certificates (PASSED as required); R4 assertion planter
(d(V)=12 in a scratch copy fires the abort).

## Independent post-run audit (assignment point 5)

`independent_audit.py` — third-path re-derivation: 11/11 pass (universe
re-enumeration with a different loop; DFS re-runs on a 44-instance sample
including all boundary classes; certificate byte-pins; both pinned checkers
replay the frozen LRAT (`s VERIFIED`, rc 0); U-witness exact expansion;
arithmetic; verdict-rule consistency; checkpoint completeness).

## Instrument defects (disclosure, pre-statement §9 + run log)

1. The first draft of the fresh CNF encoder under-constrained slot-0
   availability (missing negative unit clauses for non-input classes), which
   admitted 39/338 spurious SATs while the DFS said infeasible; found during
   pre-scouting by extracting a SAT model and replaying it exactly (slot 0 used
   a non-input operand). FIXED pre-registration (both-direction slot-0 units);
   the fixed encoder passes all direction controls and agrees with the DFS on
   all 338 instances. The broken variant was never used for any verdict.
2. drat-trim/lrat-check binaries rebuilt from the frozen source snapshot in a
   different build directory produced a different binary hash (__FILE__ path
   strings); the campaign instead freezes byte-identical rebuilt binaries whose
   SHA-256 matches the pinned record (111b0405…, b4bdebfc…).
3. Runner defect (cosmetic, fixed in place): phase decorator registered results
   after the body — KeyError on internal result bookkeeping; and two
   `subset_dfs(targets, reps)` argument mistakes were caught by the runner's
   own KeyError abort path and fixed before any verdict-affecting compute.

## exactly-one verdict

**FROZEN-NEGATIVE** for the pre-registered 55-move: the searched space (338
auxiliary directions, hash-pinned, fully enumerated, dual-instrument, zero
admitted, 73.5 s of the 3.0 CPU-h budget) is a first-class quantified negative;
sun56 is exactly 56 for its fixed orientation.
