# Pre-statement — `mws59` gap campaign: certify the MWS-59 fixed-orientation minimum (closing the widest remaining LB/UB gap)

**Created 2026-09-01, BEFORE any target compute of this campaign.**
Target: `cs/mm3/`. Gate: `mws59-gap-total59`. Agent: `Mm3NextGap` (subagent of Main).
This file is committed to git (path-scoped `git commit --only`) BEFORE
`campaign.py init`; the run dir receives a byte-identical copy with the source
commit and SHA-256 recorded.

## 0. State inherited from frozen evidence (read, not recomputed)

All machine-verified by prior frozen campaigns, cited by artifact path:

1. **Gate A/B machinery and floor instrument**
   (`campaigns/2026-08-29T235739Z_c3a32a44_7c3ed00afa49/`,
   `...010534Z_610e7116_bf045517802b/`, `...031544Z_e0f3f117_c9df97a8bf3a/`):
   the 729 Brent identities over Z for all five ladder schemes; the
   d-counter / complete memoized subset-DFS instrument
   (`gate_b_floor.py`, sha256 `c8caffa8…`); checker-certified floor UNSATs.
2. **Session-5 certified landscape**
   (`campaigns/2026-08-30T113838Z_c9532c07-61e4-40e1-8b8f-b599d1fad6dc_7e0e84591c56/certified_landscape.json`):
   `mws59|sigma^0` (all-monomial) has per-side d = (U:13, V:12, Wfac:14),
   all floors impossible, DFS states (132, 81, 672), C_lb = (14, 13, 15),
   output >= 29, total_lb = 56. The monomial-transfer theorem makes these
   values representative of the whole 48^3 monomial sandwich orbit.
3. **Session-3 artifact**
   (`campaigns/2026-08-30T012035Z_a7ea8e8e_cbdf8e63fa94/`): MWS-59 recount =
   published 59 exactly (15 left / 15 right / 29 output); a CONSTRUCTIVE
   15-addition Wfac circuit obtained by transposing their output stage is
   verified exact (`transpose_check.py`, `mws59_resolution.json`), pinning
   the upper bound C(Wfac) = 15 and output = 15 + 14 = 29.
4. **Session-8 off-diagonal census**
   (`campaigns/2026-08-31T083323Z_55447c11-5887-4fdd-aa2f-30e31f2332e8_b3046e6a6c50/`):
   all 5,796 mws59 survivor data triples decided; minimum certified total
   56, attained by exactly 2 data triples (the all-monomial one among them);
   zero data triples at <= 54.
5. **Session-9 precedent** (`20260901T123618Z_3f41e268_0264d9d39939/`): the
   identical reduction shape (stagewise exact + aux-universe existential +
   dual instrument) decided the sun56 gap. This campaign reuses that
   machinery, unchanged in semantics.

## 1. Why `mws59`, and only it, is chosen (gap table from frozen evidence only)

Certified LB vs best known for the four named families with valid fixed
orientations (frozen landscape + session-9 closure):

| decomposition/orientation | certified LB (fixed orientation) | best known | remaining gap |
|---|---|---|---|
| `paper55` σ⁰ | 55, attained (witness) | 55 | 0 — closed |
| `paper55` σ¹/σ² | 55 | 57 (circuit recomposition arithmetic) | 2, but see below |
| `sun56` σ⁰ | 56 = LB = UB (session 9) | 56 | 0 — closed |
| `sun56` σ¹/σ² | 55* (frozen shared total_lb; see §1b) | 56 by σ-covariance | 1, unprovable with this machinery |
| `mws59` σ⁰/σ¹/σ² | 56 | 59 | **3 — WIDEST** |
| `stapleton60` σ⁰/σ¹/σ² | 58 | 60 | 2, both sides open but larger state spaces |

**Chosen target: the `mws59` fixed-orientation gap [56, 59] — width 3, the
widest remaining LB/UB gap in the frozen landscape.**

### 1b. Explicit non-targets and why they are NOT decidable here (honesty list)

- `paper55` σ¹/σ² `[55, 57]`: closing the gap needs BOTH a 55-impossibility
  (one-aux existential over 3 sides — reachable, as here) AND a 56-circuit
  synthesis (two-aux scheduling along a sparsity-1 residual path — open).
  The impossibility half alone does not close the gap; deferred.
- `sun56` σ¹/σ² `[*]`: the frozen total_lb 55 for σ¹/σ² shares the σ⁰ side
  values 14/13; an LB-raising aux-1 impossibility on ONE side only lifts the
  total to 56, but the UB side needs a 56-circuit explicit synthesis in a
  σ-rotated basis. Beyond the frozen envelope's UB tools; not claimed.

  *Correction (target statement discipline): the frozen σ¹/σ² total_lb 55 is
  itself the frozen artifact's number (certified_landscape.json); this
  campaign neither strengthens nor relies on it. `[FROZEN-CITE]`*
- `stapleton60` `[58, 60]`: floor DFS state spaces 1152/1344/460 — one aux-1
  scan projected ~10-40 min/side; feasible machinery-wise, but the gap
  closes only if BOTH sides admit 59-impossibility AND a 59-circuit
  synthesis exists; width smaller than mws59's 3 → lower value for the same
  instruments.
- Non-monomial ternary reorientations (~1.0e12 orientations, four open
  directions), GL(3,Q)/GL(3,Z) beyond ternary, non-ternary alphabets,
  `laderman23_source_lock`, upper-bound synthesis for `paper55` σ¹/σ² /
  `stapleton60`: all named NOT decidable by, and excluded from, this
  machine (per frozen rule-7 boundaries).

## 2. Circuit model and addition-count semantics (fixed)

Identical to gates A/B/C (paper §2 model, same as the sun56 pre-statement §2):
- Three independent stages: left map A^9-inputs → 23 product-left-operands;
  right map B-inputs → 23 product-right-operands; output map 23 products → 9
  C-entries, counted as C(output) >= C(Wfac) + 14 via the transposition
  principle (+14 = 23 − 9, preconditions audited in-run: 23/23 nonzero
  product rows on all three sides, rank(Wfac) = 9 over Q by fraction-free
  Gaussian elimination).
- Inputs are the 9 free wires e_0..e_8 (row-major). A gate computes x+y or
  x−y of two earlier values, cost 1. Sign changes and copies cost 0.
  d(F) = number of distinct sign-classes {±t} among the target rows,
  excluding the ± input directions. Floor lemma: C(F) >= d(F); at d(F) gates
  every gate value is ± a needed class.
- The displayed upper-bound circuit of MWS-59 (their printed Table 2
  schedule: 15 left + 15 right + 29 output additions) is re-counted from the
  gate lists by direct SLP expansion (never from paper prose) in-run.

## 3. The decomposition, load direction, and data triple (fixed)

`mws59` loaded from the pinned Table-3 layout artifact
`scratch/mws59_layout.txt` lines 237–265 (three blocks; each block has 9
rows x 23 columns; column r of block 1 = the A-coefficient vector of product
r, transposed to the 23x9 U convention — same loader as the frozen
`gatec_decomps.load_mws59`). Anchors reproduced in-run before any counting:
- d(U)=13, d(V)=12, d(Wfac)=14 with frozen floor verdicts (132/81/672
  states, all floors impossible), landscape total_lb 56;
- **C(U) — BEING DECIDED, one question:** the frozen floor@13 impossibility
  gives C(U) >= 14. The registered question is the aux-1 existential at
  T = 14: does some a ∈ Au(U) make tau(U) ∪ {a} schedulable in 14 slots?
  Feasible ⇒ C(U) <= 14 ⇒ **C(U) = 14 exactly**. Infeasible ⇒ C(U) >= 15
  (exactness at 15 NOT claimed: the U-side aux-2 space is not registered).
- **C(V) — BEING DECIDED, the load-bearing side:** the frozen floor@12
  impossibility gives C(V) >= 13. Two registered existentials decide the
  ladder: aux-1 at T = 13 (tau(V) ∪ {a1} in 13 slots — feasible ⇒
  **C(V) = 13 exactly**); else aux-2 at T = 14 (pairs {a1, a2} with
  tau(V) ∪ {a1, a2} in 14 slots — feasible ⇒ **C(V) = 14 exactly**);
  else C(V) >= 15.
- **Closure for the aux-2 space (why P is the WHOLE space):** a 14-slot
  schedule for tau(V) ∪ {a1, a2} creates 14 distinct classes (12 tau + 2
  non-tau; a duplicate gate is deletable with all consumers re-pointed to
  the earlier same-class value at zero cost, so WLOG no duplicates). The
  first-created non-tau class is a signed sum over inputs ∪ tau(V) ⇒ in
  Au(V); the second (later-created) one is a signed sum over inputs ∪
  tau(V) ∪ {first} ⇒ in Au2(first). Three non-tau classes at T = 14 is
  impossible by pigeonhole (12 + 3 = 15 > 14 slots). So the registered
  existential over P = { unordered {a1, a2} : a1 ∈ Au(V), a2 ∈ Au2(a1) }
  is exactly C(V) <= 14.
- **Recursion termination:** the ladder stops at C(V) >= 15 / C(U) >= 15;
  its aux-2 (U) and aux-3 (V) continuations are NOT pre-registered — the
  verdict names the bound actually certified, never a guessed exact value.

- **C(Wfac) = 15 — ALREADY EXACT** (both directions frozen, re-verified
  in-run): d(Wfac)=14 with floor impossible ⇒ C >= 15; the
  transposed-MWS-output 15-gate circuit re-expanded exactly in-run
  ⇒ C <= 15. The upper-bound side of the gap lies in U and V, not here.


### 4b. Exact universe sizes (hash-pinned NOW, BEFORE compute)

| side | tau classes | aux-1 universe Au | csv sha256 |
|---|---|---|---|
| U | 13 | 389 | `b0d0584e616cf670c239945aa219e051aea84e4979f7cbe5084c37f05bf38f0c` |
| V | 12 | 366 | `9b0606049acf62bc595766b15d050fa049d6f1b960700268d9523fda7fbb73d1` |
| W | 14 | 421 | `5aaa81dc0a979d33f4c343fa05e2cacb2a07bf8a80af05f4671662ff8d8ed341` |

Pair universe P (V side): the SET of unordered pairs {a1, a2}, a1 ∈ Au(V),
a2 ∈ Au2(a1) := { a = canon(±x ± y) : x, y ∈ inputs ∪ tau(V) ∪ {a1} } minus
(tau ∪ inputs ∪ {a1}). Pre-registration pin: mean |Au2(a1)| ≈ 400 over a
5-sample scout, projecting |P| ≈ 79,092 — this NUMBER IS NOT LOAD-BEARING:
the scan itself derives P element-by-element over all 366 a1 (two
independent enumeration loops must agree per a1, and the per-a1 derived
set is hashed into the checkpoint), so the search space is fixed by the
construction rule, not by a projected count. The AU2 generator is closed
under canon() by construction.

Canonicalization: `canon(v) = lexicographic min(v, −v)`; inputs excluded;
zero excluded.

## 5. Instruments, arithmetic, seeds (fixed)

All arithmetic EXACT: Python ints and `fmpz` (python-flint 0.9.0); no floats
anywhere a verdict depends on them; no solver tolerance trusted for any
negative.

- **Instrument 1 (primary, complete census):** the frozen complete memoized
  subset-DFS (`gate_b_floor.py`, sha256 `c8caffa8…`) over
  classes + aux set at the given T. Complete finite search, no heuristic.
- **Instrument 2 (independent decision):** the session-9 slot-availability
  CNF encoder semantics (per-slot class, availability chaining with
  BOTH-direction slot-0 units, per-class coverage, per-slot at-most-one,
  representative (r) clauses over operand pairs), re-implemented in fresh
  code for this campaign from its section-6 spec (the sun56 encoder file
  itself is cited only for the spec; the bytes here are new), solved by
  **kissat 4.0.4** (`/opt/homebrew/bin/kissat`, one process per instance,
  `-q`). Cross-solver control instances additionally solved by
  **CaDiCaL 3.0.1** (`/opt/homebrew/bin/cadical`). The two instruments must
  agree on EVERY instance; disagreement ABORTS (defect report, no verdict).
- **Certificate layer:** every load-bearing UNSAT gets a kissat DRAT proof,
  converted to LRAT by `drat-trim -L`, accepted by BOTH pinned checkers
  (rebuilt byte-identical binaries, sha256 `111b0405…` and `b4bdebfc…`,
  matching the pinned values recorded in the session-9 campaign
  `checkers/` and `source_hashes.json`).
- Seeds: enumeration orders fully deterministic (sorted universes, fixed
  iteration order, no hash-dependent set iteration anywhere a verdict
  depends on it). Column-permutation invariance control seed 17 (frozen from
  gate C).
- Scouting disclosure: instrument-1 (DFS) runs — including FULL aux-1
  sub-space scans and sampled pair runs — were executed BEFORE this
  pre-registration to size the budget and pin §4b; instrument 2 (CNF/kissat)
  touched NO mws59 instance pre-registration. Full list and numbers in §8;
  scouting values are non-evidential and every claimed result is re-derived
  in-run under the dual-instrument protocol.

## 6. CPU / memory budget and checkpointing (fixed)

- Machine shared with sibling campaigns; all compute under `nice -n 15`,
  `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
  VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1`, `PYTHONDONTWRITEBYTECODE=1`.
- **Hard pre-registered budget: 24.0 CPU-hours** total across all
  instruments. Measured pre-registration projections (nice-15, single
  thread): V-side pair DFS 52 ms/pair ⇒ 79,092 pairs ≈ 68 min; aux-1 scans
  U 14s / V 16s / W 28s; CNF+kissat per instance ≈ 25-200 ms ⇒ pair CNF scan
  ≤ 4.4 h worst-case. Anchors, controls, witnesses, certificates ≈ 30 min.
  Projected total ≈ 6.2 h; the 24 h cap leaves 4x headroom.
- Single solver call wall cap 900 s; cumulative solver time cap 8 h; campaign
  aborts (FROZEN-INCONCLUSIVE) past either, with the exact completed prefix
  named. Memory: < 2 GB RSS (measured peak 17.2 MiB on the worst scout).
- Checkpointing: every completed instance (or a1-outer row of the pair scan)
  appends one fsynced JSON line to `run_dir/scan_checkpoint.jsonl`; the
  verdict aggregator reads only completed lines; interruption reports the
  exact completed prefix; resume-safe.
- **Budget-exhausted rule:** if the space is incomplete at cap, Verdict
  FROZEN-INCONCLUSIVE naming the completed prefix exactly; NO negative
  claim of any kind (the unfinished suffix is not evidence of absence).

## 7. Controls, both directions (fixed; all in-run AFTER init, BEFORE verdict)

**ACCEPT controls (must come out positive):**
- A1. Anchors: 729/729 Brent over Z for mws59 (pure-int AND fmpz),
  published split recount 15/15/29 = 59 by SLP expansion of the Table-2
  gate lists; certified landscape row (13/12/14 | F,F,F | 132/81/672 |
  total_lb 56) re-derived from raw factors through the same code path.
- A2. Witness ACCEPT at exact counts: the frozen transposed 15-gate
  Wfac circuit re-verified by exact expansion (covers all 23 Wfac rows);
  MWS output stage re-count 29 by direct gate-list expansion.
- A3. Extension-positive control (instrument sanity on a DIFFERENT,
  independent target — never on this campaign's own instance space):
  tau(paper55-V) — with Sun's two aux classes
  (0,0,1,0,0,1,0,0,0), (0,0,0,1,0,−1,0,0,0) — must be schedulable at
  T = 13 in both instruments (the session-9 A3 control re-run here as an
  instrument-only SAT sanity check on a known-feasible instance), and
  CaDiCaL 3.0.1 must independently confirm the kissat SAT model on that
  instance.
- A4. Frozen-row re-verification (13/12/14 d-counts, floor verdicts, state
  counts) + the transposition preconditions (23/23 nonzero rows per side,
  rank(Wfac) = 9 over Q) audited in-run.

**REJECT plants (must come out negative):**
- R1. One-addition-deleted plant: the frozen 15-gate Wfac transposed
  circuit re-serialized as an explicit SLP; with its LAST gate deleted, the
  expanded rows must mismatch Wfac (verifier REJECT by exact expansion).
- R2. Operand-perturbed plant: one operand sign flipped in the Wfac
  circuit's first gate; verifier REJECT by exact 23-row equality.
- R3. Negative-instrument plant: floor@11 on tau(V) at T = 11 must be
  INFEASIBLE in both instruments with a fresh dual-checker DRAT→LRAT
  certificate (this is a *weaker-T* instance of the frozen floor@12
  impossibility; it exercises the full certificate path in-run).
- R4. Frozen-row planter: changing one d-count assertion (d(V)=13) must
  trip the in-run abort assert, demonstrated once with an expected-failure
  run in a scratch copy, then restored.

**Both-direction guarantee:** A1–A4 AND R1–R4 are joint preconditions for
any verdict; any deviation ⇒ FROZEN-INCONCLUSIVE with the failing control
named; defects preserved and disclosed; nothing overwritten.

## 8. Scouting pre-disclosure (pre-registration honesty)

Before this file was committed, exact scouting runs WERE made, including
several FULL sub-space scans of instrument 1 (DFS) — none under instrument 2
(CNF/kissat), and no hash-pinned dual-instrument protocol was executed:
- mws59 floor DFS re-runs (d-counts and state counts reproduced: 132/81/672);
- U-side aux-1 FULL DFS scan observed: 1 of 389 feasible at T = 14;
- V-side aux-1 FULL DFS scan observed: 0 of 366 feasible at T = 13;
- W-side aux-1 FULL DFS scan observed: 8 of 421 feasible at T = 15;
- V-side pair-DFS timed on ~130 sampled pairs (0 feasible in the sample);
- universe sizes and SHA-256 pins (§4b) were computed from these scans;
- NO CNF/kissat run touched any mws59 instance pre-registration, and no
  run of either instrument ever executed under campaign.py (no init).
These scouting numbers informed the budget (§6) and the expectation
C(U) = 14, C(V) >= 14. **None of the scouting observations is evidence;
every result this campaign claims will be re-derived in-run after init
from scratch: fresh enumeration, hash-pinned universes, dual instrument,
certificates.** The scouting scripts and outputs are frozen under
PREREG_SCOUTING_NOTES/ in the run dir for the record.

## 9. Verdict rule (fixed before compute)

Let ScanU = the U-side aux-1 decision (all |Au(U)| = 389 instances at
T = 14), ScanV1 = the V-side aux-1 decision (all |Au(V)| = 366 instances at
T = 13), ScanV2 = the V-side aux-2 pair decision (all |P| = 79,092 pairs at
T = 14; |P| re-derived and asserted in-run). Every instance is decided by
BOTH instruments with per-instance agreement; any disagreement is an abort.

Branch tree (evaluated top-down; first matching branch wins):
1. **FROZEN-CERTIFIED (frontier moves DOWN, new circuit):** some ScanV1 or
   ScanV2 instance admitted by both instruments, giving an explicit
   extracted schedule; assembled with the frozen 15-gate Wfac witness and
   MWS's published left/right stage structure, the re-synthesized explicit
   total circuit has < 59 additions and passes the full 729 Brent
   identities over Z (int AND fmpz) plus per-stage exact expansion.
   The frontier move is stated as exactly what the synthesis achieves
   (e.g. "mws59 fixed orientation: UB 59 → 57"). Escalate to Main BEFORE
   any shared-ledger write. An admission whose replayed circuit does NOT
   verify (assembler defect) is NOT a negative — it is an instrument
   defect ⇒ FROZEN-INCONCLUSIVE with the defect report.
2. **FROZEN-NEGATIVE, weak branch (total >= 58):** ScanV1 = 0/366 AND
   ScanV2 = 0/|P| (both instruments, complete) AND ScanU ≥ 1 admission.
   Then C(V) >= 15, C(U) = 14 (admission ⇒ witness ⇒ exactness at 14),
   so total >= 14 + 15 + 15 + 14 = 58 and the synthesis branch failed ⇒
   mws59 fixed-orientation total ∈ [58, 59].
3. **FROZEN-NEGATIVE, strong branch (total = 59 EXACT):** ScanV1 = 0/366
   AND ScanV2 = 0/|P| AND ScanU = 0/389 (all complete, all dual-agreed).
   Then C(U) >= 15, C(V) >= 15, C(Wfac) = 15, output >= 29 ⇒
   total >= 15 + 15 + 15 + 14 = 59; MWS's published 15/15/29 = 59 witness
   (re-counted in-run) supplies total <= 59; hence
   **LB = UB = 59: Sun/MWS's published 59-addition circuit is certified
   EXACTLY optimal for its own fixed orientation**
   — the session-9-analogous closure. ScanU and the pair scan are
   load-bearing here and MUST carry certificates: every load-bearing UNSAT
   class (the floor@13/12/14 re-certs and a representative sample of ≥ 5%
   of the scanned aux instances PLUS every aux instance whose DFS state
   count exceeds 4x the frozen median) gets kissat DRAT → LRAT accepted by
   both pinned checkers.
4. **FROZEN-INCONCLUSIVE:** any control failure, any instrument
   disagreement, any certificate failure, or budget/wall exhaustion with
   the exact completed prefix named. An incomplete ScanV2 prefix admits NO
   negative inference (the unfinished suffix is not evidence of absence).
   A completed-ScanV1-negative with exhausted pair scan ⇒ INCONCLUSIVE
   naming "ScanV1 complete 0/366; ScanV2 prefix k/79,092".
- One-verdict rule: exactly one terminal verdict for the whole campaign.
  The verdict names the branch, the searched space with universes pinned
  by SHA-256, the measured CPU cost, and the exact LB/UB pair after the
  campaign.

## 10. What this campaign does NOT claim

- No claim outside the all-monomial mws59 orientation (σ⁰, (I,I,I)): the
  other 4 mws59 data triples at LB 56 and the other 10,206+ survivor data
  triples are untouched.
- No upper-bound synthesis beyond the frozen witnesses (the published
  15/15/29 circuits re-counted; no new circuit constructed unless the
  synthesis branch fires).
- No statement about `paper55`, `perminov58`, `sun56`, `stapleton60`,
  `laderman23_source_lock`, non-monomial reorientations, non-ternary
  alphabets, GL(3,ℚ)/GL(3,ℤ) beyond ternary, or other models.
- C(U) beyond the T = 14 aux-1 decision (the aux-2 U-space is not
  registered; if ScanU infeasibility alone is reached, the verdict names
  the total-58 sub-branch of §9, not 59).

## 11. Reproduction command (fixed)

```
cd /Users/jinleic/jinleic-workspace/cs/mm3/<run_dir>
OMP_NUM_THREADS=1 nice -n 15 /Users/jinleic/jinleic-workspace/cs/.venv/bin/python mws59_gap_run.py
```

with `mws59_gap_run.py` (frozen in the run dir) executing: anchors →
controls A1–A4/R1–R4 → U-side scan → V-side aux-1 scan → V-side pair scan
→ certificates → verdict aggregation, writing `verdict.json` +
`scan_checkpoint.jsonl` + `certificates/` in the run dir. Runner aborts
(exit 1) on any control failure or instrument disagreement.
