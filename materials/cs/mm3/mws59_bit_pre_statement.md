# Pre-statement — `mws59` residual bit: is the fixed-orientation minimum 58 or 59?

**Created 2026-09-02, BEFORE any compute of this campaign.**
Target: `cs/mm3/`. Gate: `mws59-bit-58vs59`. Agent: `Mm3NextGap` (subagent of
Main). Committed path-scoped BEFORE `campaign.py init`; the run dir receives a
byte-identical copy with source commit + SHA-256.

## 0. Inherited frozen facts (read, not recomputed)

From `campaigns/20260901T133118Z_8c14c3dd_7970c5318d5d/` (FROZEN-NEGATIVE,
129-file freeze, audit 7/7), for `mws59` in its fixed orientation (sigma^0,
all-monomial data triple), three-stage linear SLP model:

- **C(U) = 14 EXACTLY** — frozen floor@13 impossible (re-certified in-run with
  dual-checker DRAT→LRAT) plus an explicit 14-gate circuit extracted from the
  single admitted auxiliary of 389 and verified by exact expansion
  (`witness_U_14gates.json`).
- **C(V) >= 15** — aux-1 at T = 13: 0 of 366 admitted; aux-2 at T = 14: 0 of
  79,728 pairs admitted; both scans complete, dual instrument agreeing on every
  instance (80,483 instances total, zero disagreements).
- **C(Wfac) = 15 EXACTLY** — floor@14 impossible plus an in-run synthesized
  15-gate Wfac circuit verified by exact expansion
  (`witness_Wfac_15gates.json`).
- **C(output) >= C(Wfac) + 14 = 29** — transposition bound, preconditions
  audited (23/23 nonzero rows per side, rank(Wfac) = 9 exactly over Q).
- Therefore **total >= 14 + 15 + 29 = 58**, and the published scheme's cited
  split 15/15/29 = 59 gives total <= 59. The residual bit is 58 vs 59.
- **Defect carried forward (session-10 D1):** the frozen session-3
  `transpose_check.py` is self-circular. **No transposition-derived witness is
  trusted anywhere in this campaign**; every upper bound must be an explicit
  gate list verified by exact expansion.

## 1. The exact residual question

    total = C(U) + C(V) + C(output) = 14 + C(V) + C(output),
    C(V) >= 15 (certified),  C(output) >= 29 (certified)

**The bit is 58 ⟺ there exist BOTH a 15-addition circuit for the V-side map
AND a 29-addition circuit for the output stage.** (Then total = 58, which is
also the certified lower bound, so LB = UB = 58.) **The bit is 59 ⟺ C(V) >= 16**
(then total >= 14 + 16 + 29 = 59, and the cited published 59 closes it from
above). By the frozen scans, no other value is possible.

Equivalently in the auxiliary-ladder language of the previous campaign: a
15-gate V circuit has d(V) = 12 covered tau classes and therefore **exactly
three** auxiliary classes (fewer would strip to a 13- or 14-gate circuit, both
already refuted), so the search-theoretic form of the question is the aux-3
existential at T = 15.

## 2. Model and addition-count semantics (fixed; unchanged from the frozen gates)

Three independent stages; inputs free; a gate computes x+y or x−y of two
earlier values at cost 1; sign changes and copies free; the 23 products are not
additions; total = left + right + output. An SLP whose product operands are
inline signed sums charges (k−1) additions for a k-term operand expression and
0 for a single existing wire. Product LABELS are free: relabeling the 23
products permutes the rows of U, V and the columns of the output map
simultaneously, changes no d(F), and changes no circuit cost — this is a proved
symmetry (see §5) and is the only relabeling used.

## 3. Route P (primary, certifying): verify the published Table-2 straight-line program

`scratch/mws59_layout.txt` (the fetched arXiv:2601.05272 layout, already the
pinned source of the Table-3 factor blocks used by every frozen mws59 result)
contains the paper's **printed Table 2**: t-gates `t0..t6`, u-gates `u0..u5`,
the 23 products with inline signed operand expressions, v-gates `v0..v8`, and
the 9 output combinations. Route P transcribes it verbatim and then, IN-RUN:

1. expands every t/u gate and every product operand expression exactly over the
   9 A-entries / 9 B-entries, and every v-gate and output combination exactly
   over the 23 Table-2 products;
2. establishes the **product correspondence** to the Table-3 blocks: a
   bijection π on {0..22} and per-product signs (ε_r, ε'_r) with
   U2[r] = ε_r · U3[π(r)], V2[r] = ε'_r · V3[π(r)], ε_r·ε'_r = +1, and the
   output map matching column-wise under π with the same sign product. The
   correspondence is verified exactly; if none exists the route FAILS (no
   claim) and the campaign proceeds to Route S;
3. **recounts each stage** in the model of §2 (t/u/v gates plus inline operand
   additions plus output-combination additions);
4. verifies each stage by exact expansion: the left stage reproduces the 23
   Table-3 U rows, the right stage the 23 V rows, the output stage the 9
   C-entries as the Table-3 W combinations, all up to the free product
   relabeling π and free signs;
5. assembles the campaign circuit: **the frozen 14-gate U witness (Table-3
   orientation) + the verified Table-2 right stage + the verified Table-2
   output stage**, recounts its total, and re-verifies the underlying
   decomposition with all 729 Brent identities over Z in `fmpz` AND pure int.

**Pre-registered expectation and its falsifiable content:** the paper's cited
split is 15/15/29; the frozen session-3 machine tally agrees; a hand tally of
the printed text at pre-registration time also gives 15/15/29 (disclosed as
source reading, not compute). If the in-run recount of the right stage is 15
and of the output stage is 29, the bit is **58**. If either recount differs
from the cited value, that discrepancy is itself the finding and is reported
verbatim, with the bit decided only by what the verified gate lists support.

## 4. Route S (fallback, only if Route P fails): bounded auxiliary-3 search

**Closure argument (three auxiliary classes at T = 15).** A 15-gate circuit for
tau(V) (12 classes) creates at most 15 distinct classes; duplicates are
deletable at zero cost (rewire consumers to the earlier same-class value;
operand signs are free), and after deletion a circuit with <= 2 non-tau classes
would be a <= 14-gate circuit, refuted by the frozen scans. Hence exactly three
distinct non-tau classes a1, a2, a3 arise, and ordering them by creation time:
a1 = canon(±x ± y) over inputs ∪ tau(V) ⇒ a1 ∈ AU;
a2 over inputs ∪ tau(V) ∪ {a1} ⇒ a2 ∈ AU2(a1);
a3 over inputs ∪ tau(V) ∪ {a1, a2} ⇒ a3 ∈ AU3(a1, a2).
So the ordered-chain enumeration over (AU, AU2, AU3) is complete for the
question, and each instance is the schedulability of tau(V) ∪ {a1,a2,a3} in 15
slots.

**Pinned universe root:** AU(V) has **366** elements, csv sha256
`9b0606049acf62bc595766b15d050fa049d6f1b960700268d9523fda7fbb73d1`
(re-derived and asserted in-run; identical to the frozen campaign-1 pin).
Enumeration-only measurements taken before this commit (no feasibility
instance touched): mean |AU2(a1)| = 397.7 over a 6-element sample, mean
|AU3(a1,a2)| = 433.2 over a 24-element sample.

**Cost model (from measured frozen data).** Campaign 1's aux-2 scan decided
79,728 dual-instrument instances at T = 14 in 4,722 s = **59.2 ms/instance**.
T = 15 adds one class and one slot; the projection therefore uses 59.2 ms as a
lower estimate and 120 ms as an upper estimate per instance, and the runner
**recalibrates on its first 200 instances and records the measured rate before
any prefix boundary is declared** (measurement only; no parameter changes).

**Projected instance count:** ordered chains ≈ 366 × 397.7 × 433.2 ≈
**6.30 × 10^7**; unordered triples ≈ 1.05 × 10^7. At 59.2–120 ms that is
**173–2,100 CPU-hours**. This **exceeds the budget by two to three orders of
magnitude**, so Route S is registered as a *prefix* instrument only:

- **Tier A (certifying):** the lexicographic prefix of the ordered-chain
  enumeration in the pinned AU order (a1 ascending, then a2 ascending in
  AU2(a1), then a3 ascending in AU3(a1,a2)), dual instrument on every instance,
  per-instance fsynced checkpoint. It runs until the hard budget (§6) is
  reached and then stops at an instance boundary; the completed prefix, the
  boundary instance, the decided count, the remaining count and the measured
  rate are all frozen. **A budget-bound Tier A may claim ONLY: "no 15-gate V
  circuit exists whose auxiliary chain lies in the decided prefix". It may NOT
  claim C(V) >= 16, may NOT claim the bit, and no sampling of any kind is
  performed.**
- **Tier B (scouting, POSITIVE-ONLY, cannot certify anything negative):** a
  heuristic-ordered schedule search for a 15-gate V circuit (best-first over
  reachable classes with a node cap). A Tier-B hit is promoted only after its
  auxiliary triple is re-decided by BOTH instruments and its schedule is
  replayed as an explicit gate list and verified by exact expansion — then it
  certifies C(V) = 15 exactly like a Tier-A hit. A Tier-B miss certifies
  NOTHING and is reported as such.

## 5. Symmetry reductions — proved only

1. **canon sign-fixing** (each direction class is represented once): a gate's
   value and its negation are interchangeable because sign changes are free.
   Proved, used everywhere.
2. **Product relabeling invariance**: permuting the 23 product labels permutes
   the rows of U and V and the columns of the output map simultaneously; d(F)
   is the cardinality of a set of sign classes and is unchanged, and any
   circuit maps to a circuit of the same cost by renaming wires. Proved; used
   only to align Table-2's product numbering with Table-3's in Route P.
3. **Monomial-transfer theorem** (frozen, machine-checked in
   `gatec_invariance.py`): the verdict for this data triple transfers to its
   entire 48^3 monomial sandwich orbit. Proved (frozen); used only for the
   scope statement, never to prune a search.
4. **No other pruning is used in any certifying tier.** Tier B's heuristic
   ordering is explicitly a scouting device with no certifying power.

## 6. Budget, checkpointing, budget-bound rule

- `nice -n 15`, single thread, OMP/BLAS pinned to 1,
  `PYTHONDONTWRITEBYTECODE=1`, `sys.dont_write_bytecode`.
- **Hard pre-registered budget: 4.0 CPU-hours total.** Route P costs seconds.
  Tier B is capped at 10 CPU-minutes and 2 × 10^6 search nodes. Tier A
  receives whatever remains of the 4.0 CPU-h and **will not be extended
  mid-run under any circumstance**; the budget is enforced by a wall/CPU check
  between instances.
- Single solver call cap 900 s; memory < 2 GB.
- Checkpointing: one fsynced JSON line per decided instance (Tier A) and per
  search milestone (Tier B).
- **Budget-bound rule (restated):** the only admissible negative from a
  truncated Tier A is the prefix statement of §4; the bit remains open with
  LB 58 / UB 59 unchanged, and the verdict is a first-class MEASURED NEGATIVE
  recording instances decided, instances remaining, and the measured rate.

## 7. Controls, both directions (in-run, after init, before any verdict)

**ACCEPT:** A1 729/729 Brent over Z for the Table-3 mws59 blocks in `fmpz` AND
pure int; A2 frozen landscape row reproduced (d = 13/12/14, floors F/F/F, DFS
states 132/81/672, total_lb 56) and the frozen campaign-1 witnesses re-verified
by exact expansion (14-gate U, 15-gate Wfac); A3 extension-positive control on
an independent target (tau(sun56-V) + Sun's two aux classes at T = 13) positive
in DFS, kissat and CaDiCaL; A4 transposition preconditions re-audited.

**REJECT:** R1 one-addition-deleted plant on whichever verified stage circuit
Route P produces (or on the frozen 14-gate U witness if Route P fails) ⇒ the
exact-expansion verifier must report an uncovered target; R2 operand-sign
perturbation ⇒ verifier must reject; **R3 wrong-correspondence plant (new,
specific to Route P): a deliberately permuted product correspondence (two
labels swapped) must be REJECTED by the exact row-matching check**; R4 a
weaker-T infeasibility (tau(V) at T = 12, i.e. the frozen floor) must be
INFEASIBLE in both instruments with a fresh dual-checker DRAT→LRAT
certificate; R5 assertion planter (d(V) = 13) must abort a scratch copy.

## 8. Verdict rule (fixed before compute)

1. **FROZEN-CERTIFIED — bit = 58.** Route P (or a promoted Tier-A/Tier-B hit)
   yields a verified 15-addition V-stage circuit AND a verified 29-addition
   output stage; assembled with the frozen 14-gate U witness the total recounts
   to **58**, the decomposition passes 729/729 Brent in both arithmetics, and
   the frozen LB 58 makes it exact: **mws59 fixed-orientation minimum = 58,
   LB = UB = 58**. Consequence stated plainly: the published 59-addition scheme
   is **one addition above the optimum of its own orientation**, because its
   left stage spends 15 where 14 suffices (the campaign-1 witness); no
   published claim is contradicted — the paper claims 59 additions suffice, and
   they do.
2. **FROZEN-NEGATIVE — bit = 59.** Only if a *complete* exhaustion of the
   aux-3 space (Tier A over the whole enumeration) admits nothing: then
   C(V) >= 16, total >= 59, and with the cited/verified 59 witness the minimum
   is 59 and the published scheme is optimal for its orientation. Given §4's
   projection this outcome is not expected to be reachable inside the budget,
   and it will not be asserted on a prefix.
3. **FROZEN-NEGATIVE (measured, bounded) — bit open.** Route P fails AND
   Tier A is budget-bound: freeze the exact frontier per §6; claim only the
   prefix statement; LB 58 / UB 59 stand.
4. **FROZEN-INCONCLUSIVE** — any control failure, instrument disagreement,
   certificate failure, or a Route-P internal inconsistency that is not
   resolved into one of the above.
5. Exactly one terminal verdict; it names the route, the decided space, the
   measured cost, and the exact LB/UB pair after the campaign.

## 9. What this campaign does NOT claim

Nothing about other mws59 data triples or sigma classes, other decompositions,
non-ternary alphabets, GL(3,Q)/GL(3,Z) sandwiches, or other counting
conventions. `paper55`'s 55 remains the record and is untouched; a 58 for
mws59 is not a record and is not claimed to be one. No statement is made about
whether the paper's own toolchain would count its scheme differently.

## 10. Reproduction

```
cd /Users/jinleic/jinleic-workspace/cs/mm3/<run_dir>
OMP_NUM_THREADS=1 nice -n 15 /Users/jinleic/jinleic-workspace/cs/.venv/bin/python mws59_bit_run.py
```
Phases: A anchors + frozen-witness re-verification → B controls (A1–A4 /
R1–R5) → P Route-P transcription, correspondence, recount, assembly, Brent →
S Route-S tiers (only if P fails) → F certificates → G verdict.
