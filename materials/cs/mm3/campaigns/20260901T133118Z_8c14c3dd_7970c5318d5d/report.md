# Campaign report — gate `mws59-gap-total59` (run 20260901T133118Z_8c14c3dd_7970c5318d5d)

## Question (pre-registered: pre_statement.md, sha eb73d31f…, commit 9a90db0)

Close the certified LB/UB gap on `mws59`'s fixed orientation (Mårtensson–
Wagner–Stapleton 59-addition scheme, sigma^0 all-monomial data triple —
the widest remaining gap in the frozen landscape: certified LB 56 vs
published 59). Stagewise exact reduction: total = C(U) + C(V) + C(Wfac) + 14.

## Answer — verdict FROZEN-NEGATIVE (bounded, branch 2 of §9)

**mws59's fixed-orientation minimum moves from [56, 59] to [58, 59]:**

    total >= C(U) + C(V) + C(Wfac) + 14 = 14 + 15 + 15 + 14 = 58

with the upper bound 59 supplied by MWS's published circuit (recounted
in-run from their printed gate lists: 15 + 15 + 29 = 59). The remaining gap
is exactly one addition: 58 or 59.

Per-stage state after this campaign:

- **C(U) = 14 exactly.** Floor@13 impossible (frozen; re-certified in-run,
  `certificates/floorU_T13`, dual checkers). Aux-1 existential at T = 14
  over the hash-pinned 389-element AU(U): **1 of 389 admitted** (aux
  a = (0,0,0,−1,1,1,−1,0,1)), by BOTH instruments; the extracted 14-gate
  schedule is verified by exact expansion (`witness_U_14gates.json`:
  13 tau gates + 1 aux stepping stone, all 23 U rows covered). Witness ⇒
  C(U) <= 14; floor ⇒ C(U) >= 14; hence C(U) = 14 exactly.
- **C(V) >= 15.** Floor@12 impossible (frozen; re-certified,
  `certificates/floorV_T12`). Aux-1 at T = 13 over 366 pinned aux classes:
  **0 of 366 admitted** (dual instrument, complete). Aux-2 at T = 14 over
  the construction-rule pair universe P: **0 of 79,728 pairs admitted**
  (dual instrument, complete; P derived per-a1 by two-agreement
  enumeration, hashed per row into `scanV2_checkpoint.jsonl`).
- **C(Wfac) = 15 exactly** (both directions re-verified in-run): floor@14
  impossible (frozen; part of the frozen landscape row) + a fresh
  15-gate witness synthesized in-run from the Table-3 W block by aux-1 DFS
  (`witness_Wfac_15gates.json`, aux (0,0,−1,0,−1,1,1,0,0), verified by
  exact expansion).
- Output stage: transposition principle with preconditions audited in-run
  (23/23 nonzero product rows per side, rank(Wfac) = 9 over Q by
  fraction-free elimination) ⇒ output >= 15 + 14 = 29.

Total minimum ∈ {58, 59}; row recorded as LB 58 / UB 59. Whether 58 is
attainable requires a V-side aux-3 / U-side aux-2 search (one more spare
auxiliary slot each) — NOT registered in this campaign's space, so no
claim is made either way; and no claim is made about the other mws59 data
triples or orientations.

## Search space, budget, measured cost

- AU(U): 389 (csv sha256 b0d0584e…), AU(V): 366 (9b060604…), tau hashes
  pinned in §4b of the pre-statement and re-hashed in-run before each scan
  (mismatch ⇒ abort; none fired).
- Pair universe P: derived per a1 by the §4 construction rule; sizes
  re-derived in-run by an independent second loop on 40 probe rows (all
  agreeing); 79,728 unordered pairs decided (the pre-registration note
  projected ≈ 79,092 from a 5-sample scout — the construction rule is the
  binding definition, the projection was explicitly non-load-bearing).
- Dual instrument per instance (complete memoized subset-DFS + fresh
  slot-availability CNF solved by kissat 4.0.4), per-instance agreement
  mandatory: **agreement on all 389 + 366 + 79,728 = 80,483 instances**.
- CPU cost: 4,795.9 s ≈ 1.33 CPU-h of the pre-registered 24.0 CPU-h cap
  (single thread, nice -n 15). Peak RSS < 40 MiB. Nothing hit budget;
  the negative is COMPLETE, not budget-limited.
- Checkpointing: one fsynced JSON line per completed instance/row
  (`scanU_checkpoint.jsonl` 389 rows, `scanV1_checkpoint.jsonl` 366 rows,
  `scanV2_checkpoint.jsonl` 366 rows + per-hit rows; verdict aggregator
  and this report cite only complete files).

## Controls (both directions, all in-run post-init, all PASSED)

ACCEPT: A1 anchors 729/729 Brent over Z (int AND fmpz) + frozen rows
d/floors/states (13,12,14 | F,F,F | 132,81,672) + landscape total_lb 56;
A2 Wfac 15-gate witness re-synthesized and exact-expansion-verified; A3
extension-positive SAT control on the INDEPENDENT sun56-V target
(tau + Sun's two aux classes at T = 13: DFS True, kissat SAT, CaDiCaL 3.0.1
SAT — cross-solver agreement); A4 transposition preconditions audited.
REJECT: R1 one-gate-deleted plant REJECTED (1 tau class uncovered by exact
expansion); R2 operand-perturbed plant REJECTED (exact expansion);
R3 weaker-T tau(V) at T = 11 INFEASIBLE in both instruments with fresh
DRAT→LRAT accepted by both pinned checkers; R4 d-count planter assert
fired (rc 1, message preserved in `scratch_R4/`).

## Certificates (all dual-checker-verified by replay, rc 0 both)

- `floorU_T13` (2700 vars/… CNF sha 3f584563…), `floorV_T12` (2086 vars,
  sha pinned in verdict.json), `R3_tauV_T11`, `U_boundary_first_infeasible`
  (aux (−2,0,0,0,0,0,0,0,0) at T = 14), `V1_boundary_first_infeasible`
  (aux (−2,−2,−2,0,−2,−2,0,0,0) at T = 13). All five: `drat-trim` `s
  VERIFIED` AND `lrat-check` rc 0, binaries hash-pinned to the frozen
  2e3b2dc source snapshot (111b0405…, b4bdebfc…).

## Frontier statement (bounded, precise)

Within the three-stage linear SLP model (inputs free; gate = x±y; sign
changes free; output via the transposition bound) and for the fixed
orientation sigma^0 of the `mws59` factor blocks: **every circuit has at
least 58 additions** (C(U) >= 14 exact, C(V) >= 15, C(Wfac) = 15 exact,
output >= 29) **and the published 59-addition circuit attains 59**. The
gap [56, 59] narrows to {58, 59}. This is NOT a claim that 58 is
attainable, nor anything about other orientations, the off-diagonal
landscape, other decompositions, or other models. It does NOT contradict
any published claim: MWS's 59 stands as published; if a future 58-circuit
exists it would refute nothing frozen here (LB 58 ≤ 58).

## Instrument defects (disclosed per rule 5)

1. **D1 — the frozen session-3 transposition "constructive check" is
   self-circular.** `campaigns/2026-08-30T012035Z_…/transpose_check.py`'s
   `chain` dict is never populated (line 69 creates it; nothing appends),
   so the script's own stdout — reproduced verbatim during pre-registration
   scouting — prints `transposed gates: 0`, `Wfac circuit additions = 0`,
   `transposed circuit computes exactly the W factor map: False`, while
   `mws59_resolution.json` and the session-3 manifest claim "15 additions,
   verified to compute the W factor map exactly". The 15-addition Wfac
   upper bound was therefore UNGROUNDED by that artifact. This campaign
   re-grounds it independently: the 15-gate Wfac witness is synthesized
   from the Table-3 W block by aux-1 DFS (`witness_Wfac_15gates.json`) and
   verified by exact expansion in-run. The number 15 was CORRECT; the
   frozen verification backing it was not. The frozen campaign directory
   is untouched (append-only disclosure here; the owner may re-adjudicate
   that artifact).
2. D2 (this run, fixed pre-verdict): the primary runner's first draft
   contained a dead `transpose_output_to_wfac` path and an A3 control
   aimed at the wrong target (paper55-V with Sun's aux classes —
   infeasible by design); both were caught by phase-B smoke before any
   verdict-relevant compute, fixed (A3 re-aimed at sun56-V, the session-9
   witness-shape instance), and the smoke rerun passed. The broken drafts
   never decided any instance that enters this verdict.
3. D3 (cosmetic): a `return` statement inside the `__main__` block made
   the first hub launch crash-loop (SyntaxError at import); fixed
   pre-verdict to `sys.exit(3)` with an INCONCLUSIVE verdict file. The
   live process never reached phase C in any broken form.
4. D4 (disclosed §8 pre-registration): scouting (instrument-1-only aux-1
   scans, sample pair timings, universe hash pinning) ran before the
   prereg commit; ALL verdict-relevant numbers were re-derived in-run
   post-init under the dual-instrument protocol, and the in-run results
   match the scouting observations (U: 1/389; V1: 0/366; W-side aux-1
   8/421 was scouting context only, superseded by the in-run synthesis).
5. D5 (housekeeping): the post-run audit and witness-extraction scripts
   imported the runner via `importlib` in a shell whose
   `PYTHONDONTWRITEBYTECODE` took effect only inside the process, so a
   `__pycache__/mws59_gap_run.cpython-314.pyc` was written into the run
   dir. It contained no evidence (a byte-compiled copy of the frozen
   runner) and was deleted before `campaign.py freeze`, which refuses
   non-pristine dirs. No frozen campaign directory was touched at any
   point; the deletion is recorded here rather than silently performed.

## exactly-one verdict

**FROZEN-NEGATIVE** for the pre-registered `mws59-gap-total59` question:
the branch-2 bounded negative of pre_statement §9. mws59's fixed
orientation moves LB 56 → 58 (UB stays 59, attained by the published
circuit). The searched space is exactly the pinned AU(U) at T = 14,
AU(V) at T = 13, and the construction-rule pair universe P at T = 14 —
80,483 dual-instrument instances, zero V-side admissions, one U-side
admission whose witness was independently re-verified by exact expansion.
Measured cost 1.33 CPU-h of the 24.0 CPU-h budget. The residual open
question on this orientation — 58 vs 59 — is exactly one auxiliary slot
away on the registered recursion (V aux-3 at T = 15 or U aux-2 at T = 15
plus a synthesis) and is left to a fresh pre-registration.
