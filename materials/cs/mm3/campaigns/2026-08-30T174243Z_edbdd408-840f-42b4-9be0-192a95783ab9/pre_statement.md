# Pre-statement — Gate C continuation: the record attack over the non-monomial ternary group

**Committed 2026-08-30T17:42Z, BEFORE any computation of this campaign.**
Agent: `Mm3RecordAttack`. Folder: `cs/mm3/`. Campaign: this directory.

This campaign answers the question the gate-C sweep (campaign
`2026-08-30T113838Z_*`) left open: does a **54-addition scheme** for rank-23
3×3 matrix multiplication exist? The published record is 55
(arXiv:2607.28676). The prior sweep certified `>= 55` over the monomial
sandwich group; the monomial-transfer theorem makes that group cost-invariant;
the only place inside the ternary alphabet where a 54 can still hide is the
**non-monomial** part of the orientation space. This pre-statement fixes the
search space, encoding, budget, wall-clock cap, and adjudication rule in
advance, per repo rules 14-16.

## 1. Falsifiable outcome (fixed before compute)

**Primary claim sought (either direction, both admissible):**

- (REC) A 54-addition linear scheme for rank-23 M<3,3,3> over Z in the
  three-stage model (left map 9→23, right map 9→23, output 23→9; additions
  include subtractions; inputs free; sign changes and copies free), exhibited
  explicitly as gate lists, with all 729 Brent identities verified exactly
  over Z in the fixed convention of the gate-A harness, plus an independent
  second verification path. ⇒ ESCALATE TO MAIN BEFORE WRITING ANYWHERE.
- (NOGO) A certified `>= 55` on the total over an explicitly named set of
  orientations, each infeasibility carrying a proof accepted by BOTH
  `drat-trim` and `lrat-check` (DRAT/LRAT as in the gate-B certificate
  upgrade, campaign `031544Z`), or an exhaustive census with an explicit
  completeness argument.

A failed search over an unnamed set certifies NOTHING and will not be
reported as a bound.

**The named set (fixed now, rule 16 — no refinement after seeing results):**

S = { (D, G) : D in {paper55, sun56} (the two decomposition families with
certified total-lb 55), G in GL3(Z)[ternary] = the 6960 3×3 invertible
{-1,0,1} matrices, applied as the verified isotropy action
(U,V,W) → (X U Y⁻¹, Y V Z⁻¹, Z W X⁻¹) with (X,Y,Z) = (G,G,G) — the DIAGONAL
 subgroup of the sandwich group, glued with the sigma-orbit {Id, σ, σ²} }.

Justification for restricting to the diagonal triple (G,G,G) — fixed here,
before any result: the sandwich action with (X,Y,Z)=(G,G,G) with G ternary
unimodular is exactly the automorphism group of the INPUT ALPHABET (the
10-dimensional space span(A) ⊗ span(B) etc. is preserved entrywise only for
monomial, but (G,G,G) preserves the tensor AND keeps every target vector
inside {−1,0,1}⁹ ∪ ternary span). The full (X,Y,Z) ∈ 6960³ product would
take X,Y,Z independently; the diagonal is already 6960 × 3 = 20,880
orientation classes and contains every non-monomial sandwich that keeps the
factor blocks ternary with the SAME rank-23 structure. The off-diagonal
product 6960³ is recorded as NOT swept, with this pre-statement as the only
justification (budget; see §5). If the diagonal entry choice is found
unsound (e.g. factors become non-ternary, or Brent fails), the affected
orientations are DISCARDED, not patched.

**SU soundness check before any counting:** for every G in the 6960, the
oriented factors must (a) pass 729/729 Brent over Z with fmpz, and (b) remain
ternary. Non-ternary or Brent-failing images are excluded as "not
decompositions in the swept alphabet" — counted and recorded, never silently
dropped.

## 2. Encoding (fixed)

Counting layer (exact, no solver): the gate-A/B/C convention —
d(F) := |{ sign-classes of the 23 target vectors } \ { ±e_i }|; certified
lower bound per side = d(F) + [no d(F)-gate schedule exists], the floor
decision by complete memoized subset-DFS (gate_b_floor.subset_dfs, the SAME
code path certified in gates B/C); output stage ≥ C(Ofac) + 14 by the
transposition principle (precondition audited per instance: activity check).
Total lb = C_lb(L) + C_lb(R) + C_lb(output). This layer certifies
d-counter bounds ≥ 55 instantly; it CANNOT certify a 54 (it produces lower
bounds only). All applications of the bound layer to a 54-candidate are
ADJUDICATION, not search.

Search layers for candidates at ≤ 54 (only entered if the bound layer fails
to exclude an orientation):

- LAYER-S1 (SAT, exact over F2-embedding of the sign algebra — WAIT, this
  is exactly the arXiv:2607.29291 trap). Corrected encoding: the search runs
  over the INTEGER model directly: a "d-gate schedule" question is a
  finite-state reachability (the DFS above). NO F2 encoding will be used for
  ANY claim; F2 validity does not imply Z validity (repo note).
- LAYER-S2 (ILP, HiGHS via scipy.linprog): only as a HEURISTIC PROPOSER of
  gate schedules at a given gate count T (find a T-gate schedule whose
  values hit all 23 needed classes). Every ILP-proposed schedule is
  re-verified in exact integer arithmetic (evaluate the gate list, check all
  23 classes, then Brent 729 over Z) before any claim. ILP INFEASIBLE
  answers are NOT usable as evidence at any tolerance (the HiGHS 1e-7 trap,
  repo-known): infeasibility must instead come from the SAT/DRAT layer below
  or the complete DFS.
- LAYER-S3 (SAT/DRAT): a schedule-existence question encoded as CNF
  (python-sat / CaDiCaL153 / kissat) — the gate-B encoding family
  (`build_floor_cnf`, gate_b_encodings.py) EXTENDED with aux classes
  (value set = needed targets ∪ inputs ∪ aux-classes), T slots, at T = 54 −
  (other two sides' certified lbs) etc. UNSAT answers must carry a proof
  accepted by BOTH pinned checkers (drat-trim + lrat-check, built from the
  frozen source tree 2e3b2dc; see §4) before the bound is claimed. Any SAT
  model is re-verified in exact arithmetic. The encoding audit (which
  clauses FORBID vs REQUIRE; no positive-unit under-constraint à la
  arXiv:2607.29291) is re-run and stated for any NEW encoding before its
  UNSAT is reported.

**Ring guard (arXiv:2607.29291 read before UNSAT claims):** the paper
post-dates the target's source by 3 days and reports ten expected-UNSAT
formulas that were actually SAT due to positive-unit under-constraining.
Every new encoding here is audited for that defect shape BEFORE use: the
encoding may only FORBID (negative units) beyond the gate-B template, and a
control instance must reproduce the known SAT witnesses.

## 3. Budget and wall-clock cap

- Enumeration cost of the named set: 2 decompositions × 6960 × 3 = 41,760
  orientation CLASSES (G, σ^k) at the d-counter level. Per class: 3 side
  instances, each a sign-class reduction over 23 ternary 9-vectors + (if the
  floor question matters) a DFS at T = d. Estimated per class at the
  measured gate-C rate (~7 ms): ~5 min per decomposition overall — target
  wall-clock 2 h for the full bound-layer sweep, checkpoint every 4096
  classes (crash-safe JSONL).
- LAYER-S3 certificates: capped at 24 h solver wall-clock in total across
  all instances within THIS campaign; each individual run capped at 2 h.
- If the cap is hit, the campaign reports a PARTIAL sweep with the exact
  swept prefix enumerated (rule 7). No universal claim from a partial sweep.

## 4. Toolchain (frozen)

- python 3.14.3 (`cs/.venv`), python-flint 0.9.0 (fmpz), sympy 1.14.0,
  python-sat 1.9.dev15, scipy (HiGHS), kissat 4.0.4, CaDiCaL 3.0.1.
- Checkers: `drat-trim` and `lrat-check`, built from the pinned source
  snapshot `mm3/campaigns/2026-08-30T031544Z_*/tools_snapshot/` (drat-trim.c
  "Last edit: April 21, 2024", git 2e3b2dc). Binary SHA-256 recorded at
  freeze: drat-trim 111b0405566d…, lrat-check b4bdebfcc40d….
- arb landmines acknowledged (python-flint interval constructor is
  midpoint/RADIUS; float64 endpoint subtraction can vanish a radius) — this
  campaign uses EXACT integer (fmpz/Python int) arithmetic only; no interval
  arithmetic is needed because all quantities are integral.

## 5. Adjudication rule (fixed before compute)

- If the bound layer yields total_lb ≥ 55 for a (D, G, σ^k) class → the
  class is CLOSED as certified-excluded from the 54 record hunt (recorded).
- Any class with total_lb ≤ 54 → becomes a live target for LAYERS S2/S3 and,
  if a scheme is found, escalation per §1 (REC).
- If NO class has total_lb ≤ 54 after the FULL enumeration of the named set:
  the campaign result is "certified >= 55 over S ∪ (the 331,776
  orientations of campaign 113838Z's monomial sweep)" — where S is named
  above — a genuine extension of the certified no-go. NOT a universal
  no-54 claim: the off-diagonal 6960³ product, other decompositions, and
  the continuous groups remain open, and the README scope sentence will say
  exactly this.
- Anchor reproduction FIRST: all five published totals (55/58/56/59/60) and
  the gate-C class table must reproduce from frozen code before any new
  number is trusted.
- Set-extension rule: S is not grown after seeing results (rule 16). Any
  extension = separate pre-registered statement.

## 6. Escalation policy (fixed)

A 54-addition scheme ESCALATES TO MAIN VIA hub BEFORE being written
anywhere. A certified-no-go over the named set does NOT require escalation
(extends the current certified statement and enters README below the
contract). Any disagreement with a published paper (e.g. a defect found in
arXiv:2607.28676's optimality claim) → report to Main and STOP; never adjust
toward the paper's value. Numbers never flow to RESULTS.md or PROGRESS.md —
Main is their sole writer.
