# Pre-statement — `laderman23-ladder-total`: the exact fixed-orientation minimum of Laderman's 1976 rank-23 scheme

**Created 2026-09-02, BEFORE any feasibility or synthesis compute of this
campaign.** Target: `cs/mm3/`. Gate: `laderman23-ladder-total`. Agent:
`Mm3NextGap` (subagent of Main). Committed path-scoped BEFORE `campaign.py
init`; the run dir receives a byte-identical copy with source commit + sha256.

## 0. State inherited from frozen evidence (read, not recomputed)

- **Phase 1, `campaigns/20260902T022208Z_8c279ad1_d0276dfee63d`, verdict
  FROZEN-CERTIFIED** (git `fa3321a`): the Laderman-23 factor triple, transcribed
  verbatim from the hash-pinned primary (Bull. AMS 82(1) 1976, PDF sha256
  `a1deb200…`) and machine-verified **729/729** Brent over Z in two arithmetics.
  Triple **sha256 `522ba07f4f1784ac…`**, alphabet {-1,0,1}, 23 products.
  This campaign loads exactly that frozen JSON and re-verifies its Brent
  identities and hash in-run before any counting.
- **Gate A/B machinery and floor instrument**: `src/gate_b_floor.py`
  (d-counter + complete memoized subset-DFS), and the transposition principle
  with its preconditions, as used and frozen in sessions 5-13.
- The exact ladder for the five other named orientations (`paper55`
  sigma^0/1/2 = 55, `sun56` = 56, `mws59` = 58, `stapleton60` = 60) is frozen
  and is NOT re-opened here.

## 1. Fixed question and falsifiable outcome

**In the frozen three-stage linear-SLP model, what is the exact minimum number
of additions for Laderman's own published rank-23 decomposition at its own fixed
orientation — and how does that minimum compare to the operation count of the
scheme as Laderman printed it?**

Laderman's paper states (p.126): *"Obviously the number of additions in this
algorithm could be greatly reduced, but it is being given in its more basic
form."* So the printed count is explicitly NOT a claimed optimum. This campaign
(a) recounts the printed basic form exactly, and (b) decides the exact minimum
for that orientation, or a bounded interval.

Falsifiable both ways: a verified circuit below the certified lower bound would
refute the bound; a completed census with no admission at a level refutes the
corresponding upper-bound hypothesis.

## 2. Circuit model and addition-count semantics (fixed, identical to sessions 5-13)

Three-stage linear SLP over a commutative ring, alphabet {-1,0,1}: 9 free
inputs per side; one **gate** = `x + y` or `x - y` of two previously available
quantities and costs **1**; sign changes and copies are **free**; the 23
bilinear products are not counted. `d(F)` = number of distinct direction
classes (sign classes, `canon(v) = lexicographic min(v, -v)`) among the target
rows, excluding the +/- input directions. **Floor lemma:** every gate creates at
most one new direction class, so `C(F) >= d(F)`, with equality iff a schedule
exists in which every gate creates a needed class ("floor schedule").

**Stagewise reduction (frozen):**

    total = C(U) + C(V) + C(WFac) + 14

where `WFac[r][k]` is the coefficient of product r in output `C_k`, and the
`+14 = 23 - 9` is the output stage's transposition gap: a 9->23 circuit of L
gates transposes to a 23->9 circuit of exactly `L + 14` additions, and
conversely. Preconditions audited in-run per side: 23/23 nonzero product rows
and `rank(WFac) = 9` over Q.

## 3. Which stage counts are already exact, and which are being decided

**NOTHING is inherited as exact for this decomposition** — Phase 1 deliberately
claimed no counts. All three stage costs are being decided here, from scratch:

- `d(U) = d(V) = d(WFac) = 14` and all three floor verdicts (re-derived in-run
  by the pinned complete subset-DFS; scouted values disclosed in §8).
- Floors impossible => `C(U), C(V), C(WFac) >= 15`, hence `total >= 59`.
- **The decided question per side is whether 15 is achievable**, i.e. the
  **aux-1 existential at T = 15**: does some single auxiliary class `a` in the
  closed universe make the side's 14 needed classes plus `a` schedulable in 15
  slots? If NO for a side, that side is `>= 16`; a verified 16-gate circuit then
  makes it exactly 16.
- The output stage's attainment is decided constructively by reverse-mode
  transposition of a verified `WFac` circuit, checked entry-by-entry against
  the `WFac` map.

## 4. Closure argument bounding the auxiliary universe (the load-bearing lemma)

Fix a side with needed-class set `tau`, `|tau| = d = 14`, and suppose a circuit
uses `T = 15` gates. Each gate creates at most one new class, all 14 needed
classes must be created, so **exactly one** gate creates a non-needed class `a`
(the auxiliary), and every other gate creates a needed class. Therefore, at the
moment `a` is created, the available quantities are inputs and needed classes
only, so

    a = canon(s1 * f1 + s2 * f2),  f1, f2 in (inputs ∪ tau),  s1, s2 in {+1,-1}

i.e. **`a` lies in the signed pair-sum universe `AU(tau)`** of the 9 input
directions together with the 14 needed classes, excluding input directions and
`tau` itself. `AU` is finite and **completely enumerable**; enumerating it and
deciding each extension `tau ∪ {a}` at T = 15 therefore decides the existential
**exhaustively**, with no heuristic anywhere.

## 5. Hash-pinned universes (pinned NOW, before any feasibility compute)

Computed at pre-registration time by the closure rule of §4 (pure enumeration,
no feasibility decision), each written as sorted canonical 9-tuples, one per
line, comma-separated, hashed with sha256:

| side | d | \|AU\| | AU sha256 (first 16) |
|---|---|---|---|
| U | 14 | **408** | `a1de3df4e5077746` |
| V | 14 | **408** | `d3f3894227f939fe` |
| WFac | 14 | **408** | `13be8045ff2aec8a` |

The runner MUST re-enumerate each universe by two independent loops, check both
against these pins, and abort on any mismatch. **The searched space is
therefore exactly 3 x 408 = 1,224 dual-instrument instances**, fixed now.

## 6. Dual instruments and versions (pinned before use)

- **Instrument 1 (primary, complete census):** the frozen complete memoized
  subset-DFS `src/gate_b_floor.py::subset_dfs` (sha recorded in-run) over
  classes + auxiliary at the given T. Complete, terminating, exact; no
  heuristic, no solver tolerance.
- **Instrument 2 (independent decision):** the slot-availability CNF of the
  session-9/10 spec — vars `s[g][c]` (slot g creates class c), `a[g][k]` (base
  class k available at gate g), `r[g][c][i]` (representative choice);
  BOTH-direction slot-0 units (inputs TRUE, needed classes FALSE); per-class
  coverage; per-slot at-most-one; representative clauses over operand pairs —
  solved by **kissat 4.0.4** (`/opt/homebrew/bin/kissat`), one process per
  instance, `-q`.
- **Mandatory per-instance agreement.** The two instruments must agree on EVERY
  one of the 1,224 instances; any disagreement ABORTS the campaign with no
  verdict (FROZEN-INCONCLUSIVE) and the defect is reported.
- **Cross-solver control:** a representative instance additionally solved by
  **CaDiCaL 3.0.1** (`/opt/homebrew/bin/cadical`), which must agree.
- **Certificates:** every infeasibility claim used in a bound carries a
  kissat DRAT proof converted to LRAT by `drat-trim` and checked by **two
  independently pinned checkers** — `drat-trim` and `lrat-check`, both rebuilt
  from the frozen `2e3b2dc` source snapshot, sha256 `111b0405566d55629…` and
  `b4bdebfcc40da664be2fd416451b43f9e764e5e…` respectively (verified in-run
  before anything else runs). A certificate that either checker rejects aborts.

## 7. Symmetry reductions (fixed)

- Direction classes are canonicalized by `canon(v) = lexicographic min(v, -v)`:
  sign changes are free, so a class and its negation are one object. This is
  the only reduction applied to the search space.
- Enumeration orders are fully deterministic (sorted universes), so the census
  is reproducible independent of iteration order.
- **Column-permutation invariance control (A4 below):** permuting the 23
  products must leave `d` and every floor verdict unchanged; 8 pseudorandom
  relabelings per side, seed **17**, are checked and a mismatch aborts.

## 8. Scouting disclosure (pre-registration honesty)

Before this file was committed, exact scouting runs WERE executed on this
object, and they are why the budget below is small. Disclosed in full:

- Floor instrument on all three sides: `d = 14, 14, 14`, floor achievable
  **False** on all three, state counts 2401 each.
- Universe enumeration: `|AU| = 408` on each side (the hashes in §5 are those
  scouted values, pinned).
- Aux-1 at T = 15: **first 8 aux per side, 0 feasible**, ~50-65 ms/instance.
- Aux-2 at T = 16: 60 random pairs per side, **0 feasible**; full pair census
  measured at ~2.0 CPU-h/side (NOT pre-registered here; see §10).
- **Upper-bound scouting by a randomized greedy CSE heuristic** (120 restarts
  per side): best **16 gates on every side** (16 in 96/120 restarts). This is a
  heuristic scout; it certifies nothing.
- No CNF/SAT instance was built or solved for this object before this
  pre-registration, and no certificate was produced.

**None of the above is evidence.** Every number this campaign claims is
re-derived in-run under the dual-instrument protocol with pinned universes, and
every circuit it claims is re-verified by exact expansion independent of the
search that produced it.

## 9. Upper-bound synthesis (fixed procedure, seeds fixed)

- **Search (heuristic, allowed to be heuristic because it only proposes):**
  randomized greedy common-subexpression elimination over the side's 23 target
  rows; available set starts as the 9 free inputs; each step adds one gate
  `canon(a ± b)` from available values, staying inside {-1,0,1}; scoring =
  covered-term reduction + 3 x exact-completions + `jitter * U(0,1)`.
  **Seeds fixed now: master seed `20260902`, restart 0 with `jitter = 0`, then
  restarts 1..239 with `jitter = 2.0`, per side, in that order.** 240 restarts
  per side, deterministic given the seed.
- **Verification (the only thing that counts):** every proposed circuit is
  re-evaluated gate by gate from its operand structure over 9 symbolic free
  inputs and must reproduce all 23 target rows up to sign — an exact-expansion
  verifier that shares no code with the search. A circuit that fails
  verification is discarded and reported, never counted.
- **Output stage:** the verified `WFac` circuit is transposed by reverse-mode
  adjoint accumulation with additions counted exactly (each accumulation into a
  non-empty adjoint costs 1); the constructed stage must reproduce the `WFac`
  map in all 23 x 9 entries.
- **End-to-end:** the assembled scheme (left + right + transposed output) is
  expanded symbolically over all **81 monomials** `A_i * B_j` and must equal the
  3x3 matrix-product bilinear forms exactly, 9/9 outputs.

## 10. Budget, checkpointing, determinism

- **Hard cap 2.0 CPU-hours** total, `nice -n 15`, all BLAS/OMP thread counts
  pinned to 1, `PYTHONDONTWRITEBYTECODE=1`, `sys.dont_write_bytecode = True`
  before importing anything from a frozen campaign dir.
- Expected actual from §8 measurements: 1,224 dual-instrument instances at
  ~50-80 ms (DFS) plus ~25-200 ms (CNF) => **≈ 5-6 CPU-minutes**, plus ~30 s of
  synthesis and verification.
- **The aux-2 census is explicitly NOT in this budget.** If a side's aux-1
  census admits nothing AND no verified 16-gate circuit is found for it, that
  side's interval stays `[16, best verified]` and the aux-2 route is recorded
  as a **quantified, not-run** route at its measured ~2.0 CPU-h/side — never
  silently attempted, never claimed.
- Checkpointing: one JSONL row per instance (`side`, aux index, DFS verdict,
  CNF verdict, agreement flag, seconds), flushed as it goes, so any interruption
  leaves an exact completed prefix.
- Fully deterministic given the fixed seeds and sorted universes.

## 11. Controls — both directions, in-run, AFTER init, BEFORE any claim

**ACCEPT (must pass):**
- **A1.** The frozen Phase-1 triple re-verified in-run: hash equals
  `522ba07f4f1784ac…` and Brent **729/729** (int + fmpz).
- **A2. Known-true circuit accepted at its exact addition count.** The frozen
  `paper55` printed scheme is loaded and verified through THIS campaign's
  verifier: left 13 gates -> 23 U rows, right 14 gates -> 23 V rows, both by
  exact expansion, and the assembled total recounted as **13 + 14 + 28 = 55**.
  A known-true object must be accepted at its exact count.
- **A3.** Frozen-row re-verification: `d` and the floor verdict of one frozen
  landscape row (`paper55` U/V/WFac: d = 12/13/13, all floors impossible) must
  be reproduced by this campaign's instrument.
- **A4.** Relabeling-invariance: 8 pseudorandom permutations of the 23 products
  (seed 17) leave `d` and the floor verdict unchanged on each side.
- **A5.** Transposition preconditions per side: 23/23 nonzero rows,
  `rank(WFac) = 9` over Q.
- **A6.** Universe pins: two independent enumeration loops per side agree with
  each other and with the §5 hashes.

**REJECT (must fail):**
- **R1.** One-addition-deleted plant: delete the last gate of a verified
  16-gate circuit; the verifier must REJECT (uncovered target class).
- **R2.** Operand-perturbed plant: flip one operand sign of a consumed gate;
  the verifier must REJECT by exact expansion.
- **R3.** A floor-level infeasibility (T = 14, already certain from the floor
  DFS) re-certified through the CNF instrument as UNSAT with a **fresh
  certificate accepted by both pinned checkers**.
- **R4.** Planted wrong `d`-count assertion must fire.
- **R5.** A planted "feasible" answer: an instance known SAT (a side's classes
  at T = 15 + 5 free extra slots, trivially schedulable) must come back SAT
  from BOTH instruments — guarding against an encoder that says UNSAT always.

## 12. Verdict rule (fixed before compute) — exactly one terminal verdict

Let `ScanX` be side X's complete aux-1 census (408 instances, both instruments
agreeing per instance), and `UB_X` the smallest gate count of a
**verification-passing** circuit found for side X.

- **FROZEN-CERTIFIED** iff, for every side, the census completed and
  `LB_X = UB_X`, and the assembled scheme verified end to end. Claim: the exact
  minimum `total = C(U) + C(V) + C(WFac) + 14` for Laderman's fixed
  orientation, with the exact per-side values and the printed-basic-form
  comparison.
  - `LB_X = 15` if `ScanX` admitted and a verified 15-gate circuit exists;
    `LB_X = 16` if `ScanX` completed with zero admissions.
- **FROZEN-NEGATIVE (bounded)** iff the censuses completed but some side has
  `LB_X < UB_X`. Then the claim is exactly the interval per side and for the
  total, **naming the exact searched space (3 x 408 instances, hash-pinned) and
  the measured cost**, plus the quantified cost of the routes NOT run (aux-2 at
  ~2.0 CPU-h/side). A budget-exhausted or level-exhausted negative may claim
  **only** "no circuit of the searched form exists in the completed census" —
  never "no circuit exists", never a statement about other orientations,
  decompositions, alphabets, actions, or models.
- **FROZEN-INCONCLUSIVE** on any control failure, any instrument disagreement,
  any certificate rejection, or budget exhaustion with an incomplete census
  (naming the exact completed prefix; an incomplete prefix is NOT evidence of
  absence).
- One verdict only, for the whole campaign.

## 13. What this campaign may NEVER claim

No statement about other orientations of Laderman's tensor (only sigma^0, the
printed all-monomial orientation, is decided), other decompositions,
non-ternary alphabets, `GL(3,Q)`/`GL(3,Z)` sandwiches, other actions, other
counting conventions, or the global sub-55 question. No claim that Laderman's
printed count was wrong: the paper explicitly presents it as unoptimized, so
the comparison is a quantification of that stated slack, not a correction. No
authorship claim (Phase 1's CITED-DEPENDENCY label stands).

## 14. Reproduction command (fixed)

    cd cs/mm3/<run_dir>
    OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 15 \
      ../../../cs/.venv/bin/python laderman_ladder_run.py

Phases: A anchors (triple + Brent + frozen rows + preconditions + printed
recount) -> B controls (A2-A6, R1-R5) -> C universes (pin check) -> D aux-1
censuses (three sides, dual instrument, checkpointed) -> E synthesis +
verification -> F assembly + transposition + end-to-end -> G certificates ->
H verdict. Aborts (exit 1) on any control failure or instrument disagreement.
