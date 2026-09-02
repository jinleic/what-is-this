# Pre-statement — `sun56` gap campaign: the 55-addition question for Sun's fixed orientation

**Created 2026-09-01, BEFORE any target compute of this campaign.**
Target: `cs/mm3/`. Gate: `sun56-gap-total55`. Agent: `Mm3Sun56` (subagent of Main).
This file is committed to git (path-scoped `git commit --only`) BEFORE
`campaign.py init`; the run dir receives a byte-identical copy with the source
commit and SHA-256 recorded.

## 0. State inherited from frozen evidence (read, not recomputed)

All machine-verified by prior frozen campaigns, cited by artifact path and
SHA-256 (verified at pre-registration time):

1. **Gates A/B** (`campaigns/2026-08-29T235739Z_c3a32a44_7c3ed00afa49/`,
   `...010534Z_610e7116_bf045517802b/`, `...031544Z_e0f3f117_c9df97a8bf3a/`):
   the 729 Brent identities over Z for all five ladder schemes; the
   d-counter/floor-DFS instrument; floor-impossibility UNSATs
   checker-certified by BOTH pinned checkers (`drat-trim`, `lrat-check`,
   source snapshot `tools_snapshot/` git `2e3b2dc`, binary SHA-256
   drat-trim `111b0405…`, lrat-check `b4bdebfc…`).
2. **Gate C certified landscape** (`campaigns/2026-08-30T113838Z_…7e0e84591c56/certified_landscape.json`,
   sha256 `11fcddd1…`): `sun56|sigma^0` has per-side
   d = (U:12, V:11, Ofac:16), floors (impossible, impossible, achievable),
   C_lb = (13, 12, 16), output-stage >= 30, total_lb = 55. The monomial-transfer
   theorem (`gatec_invariance.py`, `invariance_check.json`) makes these values
   representative of the entire 48^3 monomial sandwich orbit.
3. **Off-diagonal census** (`campaigns/2026-08-30T192400Z_…0dcb5844f306/`):
   all 4,510 sun56 survivor data triples decided; min total 55, attained ONLY
   by the all-monomial triple `125_125_125`; every non-monomial triple >= 56.
   So the ONLY way a 55-total orientation of sun56 could exist inside the swept
   set is via a 55-attaining circuit on the all-monomial (i.e. unsandwiched)
   orientation.
4. **Record-attack certificates** (`campaigns/2026-08-30T174243Z_…83ab9/proofs/sun56_s0_R_d11.lrat`,
   sha256 `3c85daba…`): the floor@11 UNSAT for the sun56 right factor, accepted
   by both checkers — replayed through the pinned checkers at pre-registration
   time (`s VERIFIED` / rc 0).

## 1. The statement attacked (fixed now)

For the fixed orientation `sun56` (Sun's own factor blocks, loaded by the frozen
`sun56_verify.py` SIDES with W converted to row-major 23x9 `Wfac`), in the
three-stage linear SLP model of gates A/B/C (inputs free; gate = x+y or x−y of
earlier quantities; sign changes and copies free; total = left + right +
output additions; output counted with the transposition bound), the question is:

**Q: Does a 55-total circuit exist?**

The frozen landscape bounds every such circuit from below by
C_lb(U) + C_lb(V) + C_lb(Ofac) + 14. Two of the three sides are already exact:
- C(U) = 13: floor@12 impossible (checker-certified UNSAT, frozen) + Sun's
  13-gate U circuit, expanded exactly in-run and confirmed to reproduce all 23
  U rows over Z (witness ⇒ C(U) <= 13).
- C(Ofac) = 16: floor@16 ACHIEVABLE by the frozen DFS (16-gate explicit
  schedule, replayed in-run) ⇒ C(Ofac) = 16 exactly, and by the transposition
  principle with its preconditions audited in-run (all 23 product rows nonzero,
  exact rank-9 output map), output >= 16 + 14 = 30.
So the total minimum is exactly `13 + C(V) + 30`, and **Q ⟺ C(V) <= 12**.

**Q-decision (the entire campaign in one line):** C(V) ∈ {12, 13}.

- C(V) <= 12 ⟺ some single "auxiliary" direction a ∉ tau(V) ∪ inputs makes
  the extended class set tau(V) ∪ {a} schedulable in 12 gates (floor lemma at
  d+1: 12 gates must create 12 distinct classes = 11 tau + 1 aux; representations
  drawn from inputs and earlier-created classes only).
- C(V) >= 13 holds if no such a exists, and then Sun's 13-gate V circuit
  (witness, expanded exactly in-run) makes C(V) = 13, total = 56: **no
  55-addition circuit exists for this fixed orientation** — a first-class
  quantified negative over a precisely named finite search space.

This is a frontier move in the negative direction: `sun56` goes from **LB 55 /
UB 56** to **LB 56 = UB 56** for its fixed orientation (within the three-stage
model and the transposition accounting; NOT a universal no-55 claim about other
orientations of the tensor).

## 2. Circuit model and addition-count semantics (fixed)

Identical to gates A/B/C (paper §2 model):
- Three independent stages: left map A^9-inputs → 23 product-left-operands;
  right map B-inputs → 23 product-right-operands; output map 23 products → 9
  C-entries, bounded via the transposition principle C(output) >= C(Wfac) + 14
  where Wfac is the 9→23 "W-factor" map and +14 = 23 − 9.
- Inputs are the 9 free wires e_0..e_8 (row-major). A gate computes x+y or x−y
  of two earlier values, cost 1. Sign changes and copies are free (working with
  ± values is free). The value of a gate is a vector in Z^9 (or Z^23 for the
  output map); for the V-side question all objects live in Z^9.
- d(F) = number of distinct sign-classes {±t} among the target rows, excluding
  the ± input directions. Floor lemma: C(F) >= d(F); at d(F) gates every gate
  value is ± a needed class.
- Addition count of a displayed circuit is recomputed by direct SLP expansion
  (gate list → values), never taken from any paper text.

## 3. The decomposition and data triple (fixed)

`sun56` loaded from the frozen
`campaigns/2026-08-30T113838Z_c9532c07-61e4-40e1-8b8f-b599d1fad6dc_7e0e84591c56/scripts/sun56_verify.py`
(sha256 `3b7fb40a…`), SIDES dict; U = expand(SIDES['U']), V = expand(SIDES['V']),
W rows from expand(SIDES['W']) as 23 rows of length 9. Identity anchors
reproduced in-run before any counting:
- 729/729 Brent identities over Z with `fmpz` (exact) and independently with
  Python ints;
- 23/23 U rows, V rows reproduced from the gate lists by exact expansion;
- d(U)=12, d(V)=11, d(Wfac)=16 with frozen floor verdicts (33 / 9 / 17 DFS
  states);
- frozen landscape row values (13, 12, 16 | 30 | 55) asserted before extension.

No sandwich is applied: the data triple is the all-monomial orientation
`(X,Y,Z)=(I,I,I)`, sigma^0 — the unique lower-bound-55 data triple of the
frozen sun56 census.

## 4. Search space, complete and closed (fixed BEFORE results)

The question C(V) <= 12 is decided by exhausting a **finite, explicitly
enumerated space of auxiliary classes**:

**Aux universe Au (boxed, with the closure proof below).** Au = { canon(a) ≠ 0 :
a = ±x ± y with x,y ∈ inputs ∪ tau(V) }, canon = lexicographic min(v, −v),
excluding input classes and tau(V) members. Computed by exact integer
enumeration; size **338**, SHA-256 of the newline-joined CSV of the sorted
universe = `c10bdeeae74b601dc063cb935f7c7baca79a637cab52a7e2f85d16a46f84f813`
(repr-hash `e871c4a0…`), tau(V) CSV SHA-256 = `65e4879839d168f9…`.

**Closure proof (why 338 is the whole space, not a sample):** a 12-gate circuit
for the V-map creates 12 values; 11 distinct tau classes must appear among gate
values (each gate value is ± its class), and with 12 gates covering 11 tau
classes, at most one gate value has a class outside tau(V) — say class a
(negative values are the same class; two distinct non-tau classes would need 13
gates). Every gate value is ±x ± y of earlier values, whose classes lie in
inputs ∪ tau(V) ∪ {a}; for the aux-creating gate specifically, earlier values
have classes in inputs ∪ tau(V) (a is created once). Hence a ∈ Au by definition.
Conversely any a in Au with a schedulable extension is exactly a 12-gate circuit
(the schedule is replayed as a circuit). So **C(V) <= 12 ⟺ ∃ a ∈ Au with the
extension schedulable at 12 slots**, and the campaign computes exactly that
existential over all 338 elements. No other gate counts, no other target sets,
no sandwich, no sigma — the question is fully self-contained on the fixed V
data.

**Why this is exhaustive for Q:** section 1 reduces Q to C(V) <= 12; section 4
reduces C(V) <= 12 to the 338-instance existential. Both reductions are exact
integer logic, restated and checked in-run (assertions on gate counts and class
coverage).

## 5. Instruments, arithmetic, seeds (fixed)

All arithmetic EXACT: Python ints and `fmpz` (python-flint 0.9.0); no floats,
no interval arithmetic, no solver tolerances anywhere a verdict depends on them.

- **Instrument 1 (census):** complete memoized subset-DFS (`src/gate_b_floor.py:prep`
  + `subset_dfs`, sha256 `c8caffa8…`) applied to the extended target set
  tau(V) ∪ {a} at T = 12 for each a ∈ Au. The search enumerates every ordering
  and every representation t = ±a ± b with operands from inputs ∪ earlier
  classes; it is a complete finite census, not a heuristic. Note `prep` derives
  pair-representations from the extended class set itself, so no aux-dependent
  reachable-sum precomputation is trusted.
- **Instrument 2 (independent decision, fresh encoding):** a NEW CNF encoder
  written for this campaign (`sun56_gap_cnf.py`, frozen post-commit), s/g/c
  slot-creates-class vars, a/g/k availability vars with BOTH direction units at
  slot 0 (inputs TRUE, non-inputs FALSE — see the defect history below), chain
  propagation a[g+1][k] → a[g][k] ∨ slot-g creation, per-slot at-most-one-class,
  per-class coverage, and s → OR rep with rep → available-operand clauses. Solved
  by **kissat 4.0.4** (`/opt/homebrew/bin/kissat`, one process per instance,
  `-q`), cross-solved for control instances by **CaDiCaL 3.0.1** via python-sat
  1.9.dev15 `Cadical153`. Deterministic: no randomness beyond solver-internal
  heuristics; solver seeds pinned (`kissat` default; CaDiCaL via python-sat
  default). Both instruments must agree on every instance; disagreement ABORTS
  the campaign with a defect report.
- **Instrument 2 verification must pass both in-run direction controls
  (section 7) before any UNSAT is credited** — the pre-scout run caught the
  encoder under-constraining slot-0 availability (39 spurious SATs of 338 before
  the fix; the defect and its fix are disclosed in section 9).
- **Certificate layer:** every UNSAT credited to a verdict class also runs under
  instrument 1 (exact DFS) — the DFS is the primary certificate; additionally
  the three load-bearing UNSAT instances (floor@11 T=11; and the aux-scan
  aggregated as 338 instances at T=12; and, if reached, floor@12 for U and any
  aux-instance used in a rejected-branch summary) get explicit DRAT proofs from
  kissat, converted to LRAT by `drat-trim -L`, and accepted by BOTH pinned
  checkers from source snapshot `…/031544Z/tools_snapshot/` (git `2e3b2dc`,
  rebuilt binaries must hash `111b0405…` and `b4bdebfc…`; hashes verified at
  pre-registration on the rebuilt binaries).
- Seeds: enumeration orders are fully deterministic (sorted universes, fixed
  iteration orders). The only seed-bearing step is the column-order invariance
  control (seed 17, frozen from gate C).

## 6. CPU / memory budget and checkpointing (fixed)

- Machine shared with sibling campaigns; all compute under `nice -n 15`,
  `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
  VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1`, `PYTHONDONTWRITEBYTECODE=1`.
- Hard budget, pre-registered: **3.0 CPU-hours** total across all instruments
  (expected consumption ~0.2 h: 338 × [DFS ~30 ms + kissat ~25 ms + CNF write] ≈
  6 min, plus anchors/controls/witness replays ≈ 5 min, plus certificate
  generation for the load-bearing instances ≈ 60 min cap). Memory: < 2 GB RSS.
- Wall-clock guard: any single solver call capped at 900 s; the campaign aborts
  (FROZEN-INCONCLUSIVE) if cumulative solver time exceeds 90 min.
- Checkpointing: every completed aux index writes one JSON line to
  `run_dir/scan_checkpoint.jsonl` (crash-safe, fsync); the verdict aggregator
  reads only completed lines, so an interrupted campaign reports the exact
  completed prefix and nothing more (rule 7 compliance).
- If the budget or wall cap is exhausted before the scan completes: Verdict
  FROZEN-INCONCLUSIVE with the exact completed prefix named; NO negative claim
  of any kind (the unfinished suffix is not evidence of absence).

## 7. Controls, both directions (fixed; all run in-run AFTER init, BEFORE verdict)

**ACCEPT controls (must come out positive):**
- A1. Anchor: 729/729 Brent over Z for sun56 (fmpz AND pure-int paths), and
  published split 13/13/30 = 56 re-derived from Sun's SLP by exact expansion.
- A2. Witness ACCEPT at exact counts: Sun's V circuit verified by exact
  expansion to compute all 23 V rows; recomputed count exactly 13 (11 tau
  classes + 2 aux classes, all distinct); same for U (13) and the output stage
  (30, with Wfac reversal consistent: 30 − 14 = 16 = frozen C(Ofac)).
- A3. Extension-positive control: aux1 = (0,0,1,0,0,1,0,0,0) and
  aux2 = (0,0,0,1,0,−1,0,0,0) (Sun's own two aux classes) → the extended
  class set tau ∪ {aux1, aux2} must be schedulable at T = 13 by BOTH
  instruments (DFS True; kissat SAT), reproducing Sun's witness schedule.
- A4. Frozen-row re-verification through the same code path: the frozen
  `certified_landscape.json` sun56|sigma^0 row values (12,11,16 | floors
  F,F,T | 33,9,17 states | total_lb 55) must be reproduced in-run from raw
  factors.

**REJECT plants (must come out negative):**
- R1. One-addition-deleted plant: Sun's V circuit with one gate deleted (gate
  index 6 fixed in advance — deterministic, not cherry-picked after looking)
  must FAIL exact expansion (its final sums reference the deleted wire or
  downstream values change) or, if rewired to a still-standing 12-gate SLP,
  must fail the 729 Brent identity set; the verifier must REJECT.
- R2. Operand-perturbed plant: Sun's V circuit with one operand sign flipped
  (gate 3's second operand sign, fixed in advance) must fail either exact
  23-row equality or Brent; verifier must REJECT.
- R3. Negative-instrument plant: the floor@11 question at T = 11 must come out
  INFEASIBLE in both instruments (DFS False; kissat UNSAT with DRAT→LRAT
  accepted by both pinned checkers); and the slot-0-defect probe: the OLD
  (under-constrained) encoder variant, if resurrected for the defect appendix,
  must NOT be used for any verdict.
- R4. Frozen-row planter: changing one d-count assertion (e.g. d(V)=12) must
  trip the in-run abort assert, demonstrated once with an expected-failure run
  in a scratch copy, then restored.

**Both-direction guarantee:** A1–A4 positive AND R1–R4 negative are joint
preconditions for the verdict; any deviation ⇒ FROZEN-INCONCLUSIVE with the
failing control named, defects preserved and disclosed, nothing overwritten.

## 8. Verdict rule (fixed before compute)

- Let D = the 338-instance aux-1 decision set, instruments 1 and 2 required to
  agree per instance.
- **FROZEN-NEGATIVE for the 55-move** (the expected outcome, stated in advance
  as the prior, with immediately falsifiable content): if for ALL 338 aux
  classes both instruments decide infeasible at T = 12, then C(V) >= 13; with
  witness A2, C(V) = 13; combined with C(U) = 13, C(Ofac) = 16, output >= 30:
  every sun56 fixed-orientation circuit has **>= 56 additions**, so **no
  55-addition circuit exists for this fixed orientation within the
  three-stage linear SLP model**, and Sun's 56-addition scheme is exactly
  optimal for it. The negative is bounded exactly by: the model of section 2,
  the orientation of section 3, the aux universe of section 4 (338 elements,
  hash-pinned), and the transposition accounting. It does NOT claim: no 55
  circuit for other orientations of the same tensor; no 54 for any
  decomposition; nothing outside the fixed model.
- **FROZEN-CERTIFIED for the 55-move** (the upset outcome): if some a ∈ Au is
  admitted by BOTH instruments at T = 12, the extracted schedule is replayed as
  an explicit gate list; assembled with Sun's U circuit, the Wfac 16-gate floor
  schedule, and a 14-addition reversed-output list (from the transposition of
  Sun's own 30-addition output stage, minus the +14 gap — constructed and
  verified in-run by exact expansion), the full 55-total circuit must pass the
  729 Brent identities over Z in `fmpz` AND the independent pure-int path.
  Then the verdict is a verified 55-addition circuit for the fixed orientation,
  the frontier move is "sun56 fixed-orientation UB 56 → 55", and Main is
  messaged BEFORE anything is written to shared ledgers.
- Mixed/incomplete: FROZEN-INCONCLUSIVE (control failure, instrument
  disagreement, or budget exhaustion with the exact completed prefix named).
- The one-verdict rule: exactly one terminal verdict for the whole campaign.

## 9. Instrument defects (disclosed up front)

- D1 (found in pre-scouting, FIXED before registration): the first draft of the
  fresh CNF encoder omitted negative unit clauses forbidding non-input
  availability at slot 0; `a[0][k]` for non-input k was unconstrained, and the
  `s→OR r`, `r→a` usage chain was satisfiable by a spurious TRUE availability.
  Effect: 39/338 aux instances reported spurious SAT (e.g. aux
  (−1,0,−1,0,0,−1,0,0,1)) while the exact DFS said infeasible; extraction of a
  model and exact replay exposed slot 0 using a non-input class operand. Fix:
  unit clauses `¬a[0][k]` for every non-input k; after the fix the encoder
  agrees with the DFS on all controls and all 338 instances. The broken variant
  is preserved in this campaign for the record; no verdict rests on it.
- D2 (inherited, documented by earlier campaigns): HiGHS/ILP answers are never
  used for infeasibility claims (floating-point tolerance exposure); not used
  here at all.
- D3 (environmental): CaDiCaL153 via python-sat emits a noisy `__del__`
  traceback on interpreter shutdown (recordclass/py3.14 interaction, seen in
  campaign 031544Z); harmless, disclosed, never touches results.

## 10. What this campaign does NOT claim

No universal no-55 claim: other orientations of the same tensor, other
decompositions, and other models (e.g. different addition-count conventions,
non-ternary alphabets, basis changes) are untouched. The all-monomial
orientation's LB becoming 56 removes the last candidate for a 55-total inside
the frozen swept set of campaign `192400Z` IF combined with its census (all
other data triples >= 56); that combined statement — "the swept set contains no
55-addition scheme other than the published paper55 one, and sun56's minimum is
exactly 56" — is a COROLLARY stated only within the already-frozen census
boundary, not a new sweep.

## 11. Reproduction command (fixed)

After init/freeze, the primary artifact is re-runnable as:

```
cd /Users/jinleic/jinleic-workspace/cs/mm3/<run_dir>
OMP_NUM_THREADS=1 nice -n 15 /Users/jinleic/jinleic-workspace/cs/.venv/bin/python sun56_gap_run.py
```

with `sun56_gap_run.py` (frozen in the run dir) executing: anchors → controls
A1–A4/R1–R4 → the 338-instance dual-instrument scan → certificates → verdict
aggregation, writing `verdict.json` + `scan_checkpoint.jsonl` +
`certificates/` in the run dir. Runner aborts (exit 1) on any control failure
or instrument disagreement.
