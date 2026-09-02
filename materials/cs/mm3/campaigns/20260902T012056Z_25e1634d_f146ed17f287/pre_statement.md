# Pre-statement — `stapleton60` gap campaign: decide the fixed-orientation minimum in [58, 60]

**Created 2026-09-02, BEFORE any feasibility compute of this campaign.**
Target: `cs/mm3/`. Gate: `stapleton60-gap-total60`. Agent: `Mm3NextGap`
(subagent of Main). Committed to git path-scoped BEFORE `campaign.py init`;
the run dir receives a byte-identical copy with source commit + SHA-256.

## 0. State inherited from frozen evidence (read, not recomputed)

1. **Certified landscape** (`campaigns/2026-08-30T113838Z_…7e0e84591c56/certified_landscape.json`):
   `stapleton60|sigma^0` has per-side d = (U:14, V:14, Wfac:13), all three
   floors IMPOSSIBLE, DFS states (1152, 1344, 460), C_lb = (15, 15, 14),
   output_stage_lb = 28, **total_lb = 58**, published anchor 60.
   The monomial-transfer theorem (frozen `gatec_invariance.py`) makes these
   values representative of the entire 48^3 monomial sandwich orbit of this
   data triple.
2. **Session-8 census** (`campaigns/2026-08-31T083323Z_…b3046e6a6c50/`): all
   4,800 stapleton60 survivor data triples decided, minimum certified total
   58 (9 minimizing triples, the all-monomial one among them), zero at <= 54.
3. **Scheme source on disk**: `src/stapleton60_data.py` transcribes Stapleton's
   own Appendix-A verification script VERBATIM — `T_GATES` (6), `U_GATES` (6),
   `PRODUCTS` (23, with inline signed operand expressions), `V_GATES` (9),
   `OUTPUTS` (9). Its docstring books the published split as
   left 16 = 6 t-gates + 10 inline left additions, right 16 likewise,
   output 28 = 9 v-gates + 19 output-combination additions; total 60.
   **This campaign re-derives every one of those counts in-run by exact SLP
   expansion; no printed count is trusted.**
4. **Session-10 defect carried forward** (`campaigns/20260901T133118Z_…7970c5318d5d/report.md`
   D1): the frozen session-3 `transpose_check.py` "constructive transposition"
   is self-circular and verified nothing. **Consequence adopted here: this
   campaign never uses a transposition-derived witness. Every upper bound is
   an explicit gate list verified by exact expansion.**

## 1. Why `stapleton60` and what exactly is open

After session 10 the frozen gap table is: `paper55` sigma^0 = 55 closed;
`sun56` sigma^0 = 56 closed; `mws59` sigma^0 narrowed to {58, 59};
`paper55` sigma^1/sigma^2 [55, 57]; **`stapleton60` [58, 60]** — the target
named by the owner for this campaign. Its gap is width 2 and, unlike the
mws59 residue, it is FULLY decidable by single-auxiliary scans alone (§3).

## 2. Circuit model and addition-count semantics (fixed; identical to gates A/B/C)

Three independent stages; inputs (9 matrix entries per side) free; a gate
computes x+y or x−y of two earlier values at cost 1; sign changes and copies
free; the 23 bilinear products are not additions; total = left + right +
output additions. d(F) = number of distinct sign classes {±t} among the 23
target rows, excluding ± input directions. Floor lemma: C(F) >= d(F), and at
exactly d(F) gates every gate value is ± a needed target. Output stage is
bounded by the transposition principle C(output) >= C(Wfac) + 14 (14 = 23−9),
with preconditions audited in-run (23/23 nonzero product rows per side,
rank(Wfac) = 9 over Q by fraction-free elimination). Counting an SLP whose
product operands are inline signed sums: an operand expression with k terms
costs k−1 additions unless it is a single existing wire (cost 0); this is
exactly the model under which the frozen ladder counts 55/56/58/59/60 are
comparable, and the in-run recount uses it verbatim.

## 3. The stagewise reduction, and which counts are already exact

    total = C(U) + C(V) + C(output),    C(output) >= C(Wfac) + 14

Frozen floors give C(U) >= 15, C(V) >= 15, C(Wfac) >= 14, hence
C(output) >= 28 and total >= 58.

In-run witnesses (each an explicit gate list, verified by exact expansion —
never a printed number):
- the printed output stage expands to **28** additions and reproduces all 9
  C-entries over the 23 products ⇒ C(output) = 28 EXACTLY (>= 28 above);
- the printed left / right stages expand to **16** each and reproduce all 23
  U / V rows ⇒ C(U) <= 16, C(V) <= 16.

**Therefore total = C(U) + C(V) + 28 with C(U), C(V) ∈ {15, 16}, and the
entire remaining gap is two independent single-auxiliary questions:**

- **Q_U:** is C(U) <= 15? ⟺ ∃ a ∈ AU(U) with tau(U) ∪ {a} schedulable in
  15 slots (T = d+1 = 15).
- **Q_V:** is C(V) <= 15? ⟺ ∃ a ∈ AU(V) with tau(V) ∪ {a} schedulable in
  15 slots.

**Closure argument (why single-auxiliary is the WHOLE space at T = d+1).**
A 15-gate circuit for a 14-class target set covers all 14 tau classes, so at
most one gate value has a class outside tau (two would leave 13 gates for 14
classes). A duplicate-class gate is deletable (rewire consumers to the earlier
same-class value; operand signs are free), so WLOG the 15 gates create 15
distinct classes = 14 tau + exactly one auxiliary class a, and the gate
creating a has operands whose classes lie in inputs ∪ tau ⇒ a = canon(±x ± y)
for x, y ∈ inputs ∪ tau ⇒ a ∈ AU. Conversely any schedulable extension IS a
15-gate circuit. A zero-auxiliary 15-gate circuit would strip to a 14-gate
floor schedule, which the frozen floor decision refutes. Hence
Q_U ⟺ ∃ a ∈ AU(U) schedulable, exactly.

**Third (independent, cheap) question, pre-registered as stage 3:**
- **Q_W:** is C(Wfac) = 14 exactly? ⟺ ∃ a ∈ AU(Wfac) with
  tau(Wfac) ∪ {a} schedulable in 14 slots (T = d+1 = 14). This pins the
  Wfac map itself (not needed for the total, which already has
  C(output) = 28 exact) and replaces the void frozen transposition claim
  with a synthesized, exact-expansion-verified 14-gate witness.

## 4. Hash-pinned universes (enumerated BEFORE this file was committed; no feasibility question was touched)

| side | d | tau csv sha256 | \|AU\| | AU csv sha256 |
|---|---|---|---|---|
| U | 14 | `3b1ec60de0117cc24f2b5e7100e2deac5f72320afd39743508e3b37974528c72` | 428 | `9cb18bc76cc897e5c1de0995b37008418df062deebd034f7f4bfeeded681eb68` |
| V | 14 | `3d8092f303264f77cf4e83a81cf97a63fc0a5675c9c4b4c78d10821baa3c0dc4` | 426 | `b01988694358e18db31feb9cd3c782ca5614d1f4bd09bd6720e8fc00cf93256d` |
| Wfac | 13 | `c7c79f0aa65f0e3242d243023cd1f6ae54e42a7fcd11040c91c887a9dfa6a7c3` | 398 | `29f85ab1084e9799c0f864ca903c2c83181bc7d0abd4dca25ef092d3558f5693` |

`canon(v) = lexicographic min(v, −v)`; zero and ± input directions excluded.
Universes are re-derived in-run and the hashes asserted before each scan
(mismatch ⇒ abort). **Enumeration-only pre-commit disclosure:** the three
universes above were enumerated (deterministic set arithmetic over the frozen
factor blocks) to pin these hashes; **no DFS decision, no CNF, and no
feasibility instance of any kind was run before this commit** — unlike the
session-10 campaign, whose full aux-1 scans pre-dated its prereg and were
disclosed as such.

## 5. Instruments, arithmetic, symmetry, seeds

- **Instrument 1 (primary, complete census):** frozen complete memoized
  subset-DFS `src/gate_b_floor.py` (sha256 `c8caffa8…`).
- **Instrument 2 (independent):** the slot-availability CNF of the
  session-10 runner (per-slot class vars, availability chain, BOTH-direction
  slot-0 units, per-class coverage, per-slot at-most-one, representative
  clauses over operand pairs) solved by **kissat 4.0.4**; control instances
  cross-solved by **CaDiCaL 3.0.1**. Per-instance agreement MANDATORY on
  every instance of every scan; any disagreement aborts the campaign.
- **Certificates:** every load-bearing UNSAT gets kissat DRAT → `drat-trim -L`
  LRAT → accepted by BOTH pinned checkers (`drat-trim` `111b0405…`,
  `lrat-check` `b4bdebfc…`, from the frozen 2e3b2dc snapshot). Load-bearing
  set: the three floor instances (U@14, V@14, Wfac@13) and, for each scan
  that returns all-infeasible, the lexicographically first infeasible
  instance plus a seeded sample of 10 (seed 20260902).
- All arithmetic exact (Python int and `fmpz`); no float, no tolerance.
- **Symmetry reductions:** (i) canon() sign-fixing (each direction class
  appears once); (ii) the frozen monomial-transfer theorem — the verdict for
  this data triple transfers to its whole 48^3 monomial sandwich orbit, so no
  sandwich is enumerated; (iii) no other reduction is applied and none is
  needed: the scans are complete over the pinned universes.
- Seeds: enumeration orders deterministic (sorted); the only seeded step is
  the certificate sample (seed 20260902).

## 6. Budget, checkpointing, and the budget-bound rule

- All compute `nice -n 15`, single thread, `OMP_NUM_THREADS=1` (and the four
  sibling BLAS vars), `PYTHONDONTWRITEBYTECODE=1`, `sys.dont_write_bytecode`.
- **Hard pre-registered budget: 6.0 CPU-hours.** Measured session-10
  analogues project: aux-1 instance ≈ 30–200 ms dual-instrument, so
  428 + 426 + 398 = 1,252 instances ≈ 1–5 min; anchors, recounts, controls,
  certificates, assembly and Brent ≈ 5 min. Expected total < 15 min; the cap
  is 24x headroom. Single solver call cap 900 s; cumulative solver cap 2 h.
  Memory < 2 GB.
- Checkpointing: one fsynced JSON line per decided instance into
  `scan_<side>_checkpoint.jsonl`; the aggregator reads only complete lines.
- **Budget-bound negative rule:** if a scan does not complete, the campaign
  may claim ONLY (a) the frozen floor lower bounds, (b) the in-run verified
  printed-witness upper bounds, and (c) that the completed lexicographic
  prefix of the named universe contains no admissible auxiliary — never that
  the side value is 16, and never a total. Verdict then
  FROZEN-INCONCLUSIVE with the exact prefix length and measured cost.

## 7. Controls, both directions (all in-run, after init, before any verdict)

**ACCEPT (must be positive):**
- A1. 729/729 Brent identities over Z for stapleton60 in `fmpz` AND pure int;
  ternary check; frozen landscape row reproduced through the same code path
  (d = 14/14/13, floors F/F/F, DFS states 1152/1344/460, total_lb 58).
- A2. Printed-SLP recount by exact expansion: left 16, right 16, output 28,
  total 60, each stage reproducing all its target rows exactly (23 U rows,
  23 V rows, 9 C-entries over the 23 products).
- A3. Extension-positive control on an INDEPENDENT target: tau(sun56-V) with
  Sun's two auxiliary classes at T = 13 must be schedulable in BOTH
  instruments, and CaDiCaL must confirm the kissat SAT.
- A4. Transposition preconditions audited: 23/23 nonzero rows on all three
  sides, rank(Wfac) = 9 exactly over Q.

**REJECT (must be negative):**
- R1. One-addition-deleted plant: delete the last tau-creating gate of the
  in-run Wfac (or printed-left) witness ⇒ the exact-expansion verifier must
  report an uncovered target class.
- R2. Operand-perturbed plant: flip one operand sign of a consumed gate ⇒
  the verifier must reject by exact recomputation.
- R3. Negative-instrument plant: a weaker-T instance (tau(Wfac) at T = 12)
  must be INFEASIBLE in both instruments with a fresh dual-checker
  DRAT→LRAT certificate.
- R4. Assertion planter: a wrong d-count assertion (d(U) = 15) must abort a
  scratch copy.

Any control deviation ⇒ FROZEN-INCONCLUSIVE naming the failing control;
defects preserved and disclosed.

## 8. Verdict rule (fixed before compute)

Let ScanU (428 instances at T = 15), ScanV (426 at T = 15), ScanW (398 at
T = 14), all dual-instrument with mandatory per-instance agreement.

1. **FROZEN-CERTIFIED — gap closed downward.** If ScanU and/or ScanV admit,
   the admitted schedules are replayed as explicit 15-gate circuits, verified
   by exact expansion independent of the search, and assembled with the
   verified printed stages into a full scheme whose total is
   C(U) + C(V) + 28 ∈ {58, 59}; the assembled scheme must pass all 729 Brent
   identities over Z in `fmpz` AND pure int. Claim: stapleton60's fixed
   orientation has minimum EXACTLY that total (LB = UB), improving the
   published 60 by 2 or 1. Escalate to Main before any shared-ledger write.
2. **FROZEN-NEGATIVE — published optimum certified.** If ScanU and ScanV are
   both complete with zero admissions, then C(U) = C(V) = 16 (floors + the
   in-run verified 16-gate printed witnesses) and total = 16 + 16 + 28 = 60
   EXACTLY: **Stapleton's published 60-addition scheme is optimal for its own
   fixed orientation**, and the frozen LB 58 was not tight. Bounded exactly
   by the model (§2), this orientation and data triple (§0.1), and the pinned
   universes (§4).
3. **Mixed** (one side admits, the other is complete-negative) is covered by
   branch 1 with the exact total 59.
4. **FROZEN-INCONCLUSIVE** — any control failure, instrument disagreement,
   certificate failure, or budget-bound incompleteness (§6 rule).
5. ScanW is reported in every branch: admission ⇒ C(Wfac) = 14 exactly with
   a verified 14-gate witness; all-infeasible ⇒ C(Wfac) >= 15, which would
   RAISE output to >= 29 and hence total to >= 59 — a tension with the
   in-run recount output = 28 that would be escalated as an instrument
   defect, not silently absorbed.
6. Exactly one terminal verdict; it names the branch, the pinned searched
   space, the measured CPU cost, and the exact LB/UB pair after the campaign.

## 9. What this campaign does NOT claim

Nothing about other stapleton60 data triples or orientations (the census
covers their LBs only), other decompositions, non-ternary alphabets,
GL(3,Q)/GL(3,Z) sandwiches, other addition-count conventions, or the
`mws59` 58-vs-59 residue. No universal statement about rank-23 3x3
additive complexity; `paper55`'s 55 remains the record and is untouched.

## 10. Reproduction

```
cd /Users/jinleic/jinleic-workspace/cs/mm3/<run_dir>
OMP_NUM_THREADS=1 nice -n 15 /Users/jinleic/jinleic-workspace/cs/.venv/bin/python stapleton60_gap_run.py
```
Runner phases: A anchors/recount → B controls → C ScanU → D ScanV →
E ScanW → F certificates → G assembly+Brent → H verdict. Aborts (exit 1) on
any control failure or instrument disagreement.
