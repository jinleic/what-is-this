# Pre-statement — `paper55` sigma^1 / sigma^2 orientations: decide each minimum in [55, 57]

**Created 2026-09-02, BEFORE any compute of this campaign.**
Target: `cs/mm3/`. Gate: `paper55-sigma12-total55`. Agent: `Mm3NextGap`
(subagent of Main). Committed path-scoped BEFORE `campaign.py init`; the run dir
receives a byte-identical copy with source commit + SHA-256.

## 0. Frozen inputs (read, not recomputed)

1. **Certified landscape** (`campaigns/2026-08-30T113838Z_…7e0e84591c56/certified_landscape.json`),
   for `paper55`:
   | class | per-side d | floors | C_lb | output_lb | total_lb | best known |
   |---|---|---|---|---|---|---|
   | sigma^0 | (12, 13, 13) | F, F, F | (13, 14, 14) | 28 | **55** | 55 (attained) |
   | sigma^1 | (13, 13, 12) | F, F, F | (14, 14, 13) | 27 | **55** | 57 |
   | sigma^2 | (13, 12, 13) | F, F, F | (14, 13, 14) | 28 | **55** | 57 |
   DFS state counts (116, 66, 33) for sigma^1 and (66, 33, 116) for sigma^2.
2. **Gate B (frozen, checker-certified)**: for `paper55` sigma^0 all three side
   costs are EXACT — **C(U) = 13, C(V) = 14, C(Wfac) = 14** — each by
   floor-impossibility at d(F) (DRAT/LRAT accepted by both pinned checkers in
   campaign `…031544Z…`) plus a witness at d(F)+1.
3. **Gate A (frozen)**: the paper's printed straight-line program is on disk as
   explicit gate lists — `src/tensor_data.py` `LEFT_SLP` (13 gates),
   `RIGHT_SLP` (14 gates), `OUTPUT_SLP` (28 gates) with `C_ALIASES` — verified
   729/729 Brent over Z with counts 13/14/28 = 55.
4. **Frozen sigma action** (`gatec_sweep.sigma_orbit`): with `T` = per-product
   3x3 transpose of the 9-vector rows,
   `sigma^1 = (V, T(W), T(U))`, `sigma^2 = (T(W), U, T(V))`.
   All 15 landscape classes are Brent-verified 729/729 (`brent_failures: 0`).
5. **Defect carried forward (session-10 D1)**: no transposition-derived
   *claim* is ever evidence here. Where a transposition CONSTRUCTION is used,
   the evidence is the resulting explicit gate list plus an end-to-end exact
   expansion — never the principle.

## 1. Exact orientation definitions relative to the sigma^0 record

Writing (U, V, W) for the frozen sigma^0 blocks (23x9 row-major, W in Wfac
rows), the two target orientations are exactly the frozen ones:

- **sigma^1** = (V, T(W), T(U)): left map = V, right map = T(W), output-factor
  map = T(U);
- **sigma^2** = (T(W), U, T(V)): left map = T(W), right map = U,
  output-factor map = T(V).

`T` is a fixed permutation of the 9 input coordinates (the involution
`[0,3,6,1,4,7,2,5,8]`, i.e. transpose of each 3x3 block).

## 2. Model and addition-count semantics (unchanged, fixed)

Three stages; inputs free; gate = x+y or x−y at cost 1; sign changes and copies
free; 23 products not counted; total = left + right + output additions; an
SLP whose operands are inline signed sums charges (k−1) for a k-term operand.
`C(output) >= C(Ofac) + 14` (transposition bound, 14 = 23−9), preconditions
audited in-run (23/23 nonzero rows per side, rank of the output factor 9 over Q).

## 3. Stagewise reduction: which counts are already exact per orientation

**Proved symmetry (used here, stated as a theorem to be re-checked in-run):**
applying a fixed permutation/sign pattern to the input coordinates carries any
circuit to a circuit of the same cost (rename the free input wires), and
permuting the 23 product labels likewise costs nothing. Hence for any block X,
`C(T(X)) = C(X)`, and a circuit for X becomes a circuit for T(X) by input
renaming alone.

Therefore, from the frozen gate-B exact values:

| orientation | left | right | output-factor | resulting bound |
|---|---|---|---|---|
| sigma^1 | C(V) = **14** exact | C(T(W)) = C(Wfac) = **14** exact | C(T(U)) = C(U) = **13** exact | output >= 27, total >= 14+14+27 = **55** |
| sigma^2 | C(T(W)) = **14** exact | C(U) = **13** exact | C(T(V)) = C(V) = **14** exact | output >= 28, total >= 14+13+28 = **55** |

This reproduces the frozen total_lb 55 for both classes, now with **every side
cost exact rather than a floor+1 bound**. Consequently the ONLY open quantity
per orientation is whether its output stage attains its bound (27 for sigma^1,
28 for sigma^2). The campaign therefore decides each orientation by
constructing and verifying the missing explicit circuits.

## 4. What will be constructed and how it is verified (no search unless step 3 fails)

For each orientation the campaign assembles three explicit stage circuits and
verifies each by exact expansion, then verifies the whole scheme end to end:

- **left / right stages**: reuse the frozen printed circuits under input
  renaming — `RIGHT_SLP` (14 gates) for sigma^1's left map; `LEFT_SLP` (13
  gates) for sigma^2's right map; and for the `T(W)` side a **14-gate Wfac
  circuit synthesized in-run** by the aux-1 instrument at T = d(Wfac)+1 = 14
  over the hash-pinned auxiliary universe (the frozen record states aux-1 at
  d+1 is possible for all three paper55 blocks; this campaign re-derives it and
  keeps the explicit gate list). **No transposition-derived Wfac witness is
  used.**
- **output stage**: constructed by **reverse-mode transposition** of the
  corresponding verified output-factor circuit (the 13-gate `LEFT_SLP` for
  sigma^1, the 14-gate `RIGHT_SLP` for sigma^2, each under input renaming),
  implemented here as an explicit gate list. The construction is NOT evidence;
  the gate list is counted directly and verified by exact expansion to compute
  the 9 outputs as the orientation's W-combinations of the 23 products.
- **end-to-end**: the assembled scheme is expanded symbolically over all 81
  monomials A_i·B_j and compared, entry by entry, against the bilinear map
  defined by that orientation's frozen factor triple (equivalently the matrix
  product in that orientation's convention), and the orientation's factor
  triple is re-verified with all 729 Brent identities over Z in `fmpz` AND
  pure int.

Expected counts: sigma^1 = 14 + 14 + 27 = **55**; sigma^2 = 14 + 13 + 28 =
**55**. These are *predictions*, not assumptions: the recount comes from the
explicit gate lists, and if a constructed output stage costs more than its
bound the campaign reports the measured count and the resulting UB instead.

## 5. Search component, closure argument and hash-pinned universe

The only search is the Wfac aux-1 synthesis (needed twice, same instance set):
`C(Wfac) <= 14` ⟺ some a in AU(Wfac) makes tau(Wfac) ∪ {a} schedulable in 14
slots. Closure: a 14-gate circuit for a 13-class target set has at most one
non-tau gate value (two would leave 12 gates for 13 classes), duplicate-class
gates are deletable at zero cost (operand signs free), and a zero-auxiliary
14-gate circuit would strip to the frozen-refuted 13-gate floor schedule; the
auxiliary's creating gate has operands with classes in inputs ∪ tau, so
a ∈ AU. Universe pinned by construction rule and by hash in-run:
**|AU(Wfac_paper55)| and its csv sha256 are asserted in-run against the value
derived by two independent enumeration loops** (the pre-registration does not
hard-code a number it has not enumerated; the pin is the construction rule plus
the two-loop agreement, and both are recorded in the checkpoint).

Instruments: the frozen complete memoized subset-DFS (`src/gate_b_floor.py`,
sha256 `c8caffa8…`) and the slot-availability CNF solved by **kissat 4.0.4**,
with **per-instance agreement mandatory**; control instances cross-solved by
**CaDiCaL 3.0.1**; every load-bearing UNSAT gets kissat DRAT → `drat-trim -L`
LRAT accepted by BOTH pinned checkers (`111b0405…`, `b4bdebfc…`).

## 6. Proved symmetry reductions (only these)

1. canon sign-fixing (free negation).
2. **Fixed-coordinate relabeling**: input permutations/sign patterns are free ⇒
   `C(T(X)) = C(X)`; re-checked in-run by recomputing d and the floor verdict
   for T(X) and comparing with X.
3. **Product relabeling** is free (permutes rows of all three blocks
   simultaneously).
4. Frozen **monomial-transfer theorem** for orbit transfer of the verdict.
No other pruning; nothing heuristic enters a certifying step.

## 7. Cost model, projected counts, hard budget

Measured rates from this session's frozen campaigns: aux-1 dual-instrument
instances cost 30–200 ms (stapleton60: 1,252 instances in 110.7 s = 88 ms;
mws59 T=14 pairs: 59.2 ms). The Wfac aux-1 universe for paper55 is expected in
the 350–450 range by analogy with the frozen sides (398–428 measured on
stapleton60, 366–421 on mws59) ⇒ projected **< 2 CPU-minutes**. All other work
(expansion, transposition construction, recounts, 729 Brent x 2 arithmetics x 2
orientations, controls) is arithmetic: seconds.

**Hard pre-registered budget: 2.0 CPU-hours**, never extended mid-run;
`nice -n 15`, single thread, OMP/BLAS pinned to 1, `PYTHONDONTWRITEBYTECODE=1`.
Per-instance fsynced checkpointing for the scan. No sampling anywhere.

## 8. Verdict rule (fixed before compute)

Per orientation (sigma^1, sigma^2), independently:

1. **CERTIFIED-55**: all three stage circuits verified by exact expansion, the
   recount totals **55**, and the end-to-end symbolic expansion matches that
   orientation's bilinear map on all 81 monomials with 729/729 Brent in both
   arithmetics. Then that orientation's minimum is **exactly 55** (LB 55 frozen
   and re-derived with exact side costs; UB 55 attained), superseding its
   "best known 57" row. Stated plainly: this is the **record value at a
   different orientation of the same tensor**, not a new record.
2. **CERTIFIED-ABOVE-55**: if a constructed output stage verifies but costs
   more than its bound, the campaign reports the exact verified total T and the
   exact minimum interval [55, T]; if in addition the frozen LB is shown
   attained-or-not by a completed search it states the exact minimum. Any claim
   that an orientation is strictly worse than the record orientation must name
   the exact certified minimum, not an interval.
3. **BUDGET-BOUND / INCOMPLETE**: if the Wfac aux-1 scan does not complete, the
   only admissible statement is the decided prefix; no side value, no total,
   and the orientation's interval stays [55, 57]. Verdict FROZEN-INCONCLUSIVE
   with instances decided, remaining, and measured rate.
4. **FROZEN-INCONCLUSIVE** also on any control failure, instrument
   disagreement, certificate failure, or an internal inconsistency.
5. One terminal verdict for the campaign covering both orientations; it names
   each orientation's outcome, the searched space, measured cost, and the exact
   LB/UB pair per orientation.

## 9. Controls, both directions (in-run, after init, before any verdict)

**ACCEPT** — A1: 729/729 Brent over Z (`fmpz` AND pure int) for sigma^0,
sigma^1 and sigma^2 triples; A2: the frozen printed circuits re-verified by
exact expansion (LEFT_SLP → 23 U rows at 13 additions; RIGHT_SLP → 23 V rows at
14; OUTPUT_SLP + C_ALIASES → the 9 outputs at 28) and the sigma^0 total
recounted as 55; A3: the frozen landscape rows for sigma^1/sigma^2 reproduced
through the same code path (d triples, floor verdicts, DFS state counts
116/66/33 and 66/33/116); A4: the proved relabeling symmetry re-checked
numerically (d and floor verdict of T(X) equal those of X for X in {U, V, W});
A5: transposition preconditions audited per orientation.

**REJECT** — R1: one-addition-deleted plant on each verified stage circuit ⇒
exact-expansion verifier must reject; R2: operand-sign perturbation ⇒ must
reject; R3: a wrong-orientation plant — assembling sigma^1's stages against
sigma^2's factor triple must FAIL the end-to-end expansion; R4: the frozen
floor instance for Wfac at T = 13 must be INFEASIBLE in both instruments with a
fresh dual-checker certificate; R5: assertion planter (d(Wfac) = 14) must abort
a scratch copy.

## 10. What this campaign does NOT claim

No new record: `paper55` sigma^0's 55 stands as the record and a 55 at
sigma^1/sigma^2 is the same value at a different orientation of the same
tensor. Nothing about other decompositions, other data triples/sandwiches,
non-ternary alphabets, GL(3,Q)/GL(3,Z), other counting conventions, or the
global sub-55 question.

## 11. Reproduction

```
cd /Users/jinleic/jinleic-workspace/cs/mm3/<run_dir>
OMP_NUM_THREADS=1 nice -n 15 /Users/jinleic/jinleic-workspace/cs/.venv/bin/python paper55_sigma_run.py
```
Phases: A anchors + frozen-row + printed-circuit re-verification → B controls
(A1–A5 / R1–R5) → C Wfac aux-1 synthesis (dual instrument, checkpointed) →
D sigma^1 assembly + end-to-end expansion → E sigma^2 assembly + end-to-end
expansion → F certificates → G verdict.
