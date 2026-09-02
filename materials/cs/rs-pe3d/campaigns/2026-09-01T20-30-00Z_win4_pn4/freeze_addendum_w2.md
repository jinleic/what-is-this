# Freeze addendum — campaign `2026-09-01T20-30-00Z_win4_pn4` (Campaign W, terminal)

Appended below the existing gate contract (append-only). UTC at writing:
2026-09-01T05:2xZ. Agent: RsPe3dW2 (finishing agent for the RsPe3dMW handback).
Grandfathered-legacy campaign per the 2026-09-01 workspace control plane:
no migration, no renaming, frozen bytes untouched — everything below is
additive files plus this append.

## A. Result (frozen, this campaign's own instruments)

**Census (PARITY PRIMARY, pre-registered showcase):** the weight-EXACTLY-4
window at PN4 = (421,(4,5,7),t=(1,1,1),Λ=Id) was censused COMPLETELY by the
parity instrument: all C(140,4) = 15,329,615 supports, exact F_421
arithmetic, emptiness test dim(V∩F^S) = |S| − rank(π[:,S]), π = H0⊗H1⊗H2
(72×140), column index (a0·s1+a1)·s2+a2 — crucially including the
non-candidate side (a support is settled as a non-candidate by rank = 4,
with no k-per-subset test), so no missed-candidate failure mode exists by
construction.
- **Candidates: EXACTLY 35** (`parity_w4_ledger.jsonl` kind "final",
  next_index 15329615, hits_total 35; full per-hit detail in
  `parity_w4_hits.json`).
- **The candidate set is SET-EQUAL to the 35 direction-0 minimum lines**
  {o, o+35, o+70, o+105}, o ∈ [0,34] (`parity_w4_hits.json` S fields vs
  src.rs Inst line_indices; 0 non-line supports, 0 missing lines, 0
  duplicated candidates). Killing the owner-input status: this is OUR
  number, machine-verified, from our own kernel.
- All 35 have **dim(V∩F^S) = 1**; each candidate class is the single
  normalized all-ones word on its line.
- Engagement-rate figure (in-process, this run): first 2.4M supports in
  ~81 s ≈ 31,500 supports/s (supports/s; instrument = accounting does not
  enter the mathematical claims).
- Direction-1/2 minimum lines (28 of size 5, 20 of size 7): each is
  dependent as a FULL line with dim 1 and every 4-subset independent
  (rank-4 of a 5/7-set; `parity_line_checks()` printed in
  `parity_run.log` and part of `parity_w4_ledger.jsonl`'s processes... see
  the addendum file `w4_owner_input_checks.json` for the rerun record).
  Consistent with d1=5, d2=7: consistent with the owner INPUT, not
  re-labeled as our own independent full-line sweep (we did reproduce the
  rank facts).

## B. Set comparison — what is element-wise and what is honestly partial

Per the handback contract, the V-basis ladder (`w_ladder`, pid 85270,
census_runner.py --anchors-only) was left running to completion. At
freeze-writing time it had landed **B01-B10 exact** (55/20/55/110/30/110/
24/280/16/40, tally-matched 13-gate-B rows plus this one) and was **still
inside B11** (~2.3M-supports row > 1h). B11/B12/B13 + the B12/B13 blocks
had NOT landed at freeze time; census_ledger.jsonl remains at the kill-demo
checkpoint (2 lines, idx 410000).

What IS element-wise on disk and frozen:
- `w4_vbasis_verification.json`: the **V-basis kernel re-checked all 35
  parity-PRIMARY candidates** — every one nonzero (rank(B[:,off]) < 68),
  every dim = 1 matching the parity dim. **0/35 disagreements** on the
  found set.
- The V-basis kernel's full-space sweep over all 15,329,615 (which alone
  could surface a parity-INVISIBLE extra candidate) did NOT complete this
  session (two handback-era background copies were killed by the handback
  itself; frozen predecessor figure 3,449 s; per the owner's 465/s
  ladder rate ≈ 80 min/row). The V-basis FULL census is therefore
  recorded as **NOT COMPLETED in this session** — honestly disclosed;
  the ladder process keeps running and its rows land into
  `census_ledger.jsonl` post-freeze (append-only; a follow-up freeze can
  fold them in).
- Per the contract, "two kernels both reporting 35" was not achieved
  end-to-end; what was achieved is: parity PRIMARY complete (including
  the non-candidate side by construction), V-basis nonzero-check on the
  full 35-candidate set (agreement 35/35, dims matched), and the V-basis
  kernel pre-anchored on 10/13 gate-B rows + P2 anchors (55/20/55/110/30/
  110/24/280/16/40 = tallies, PN6 3486, PP2 917) re-landed exactly
  (anchors_run.log; parity ledger cross-check 15/15 MATCH).

## C. Rho stage (step 3 — exact closure over the candidate set)

`rho_w4_results.json` (frozen). All 35 candidates: **delta(S) = 4
exactly**, closed by the frozen exhaustive ascending-cost sweep
(DeltaEngine.delta_exact, unchanged src) with line_tuple = a SINGLE
direction-0 line in every case (the membership test: the all-ones word on
one d0-line lies in W(B) for B = that line). **ρ^window_w4 = min wt/delta
= Fraction(4,4) = 1/1 = 1** exactly. No candidate has delta > 4 (none
needed a costlier closure) and none has delta < 4 (wt = 4 lower-bounds
delta? no — delta < 4 is the sweep's own business: none was found).
Per-candidate rho values: all 35 equal 1 (field "rho" = "1", i.e. "4/4").

## D. Verdict (step 4)

- **Pattern P at PN4, w=4 window:** P makes NO committed prediction here
  (its pre-registered scope is the ≤ 3-window; e.g. P2's arm rows all
  have min(s) ≤ 3). What we OBSERVE is ρ^window_w4 = 1 — the "all
  line words" floor outcome flagged as expected in the pre-statement
  (rho = 1 branch, "consistent with the ≤3-window history at non-{2,3}
  rows"). Retrospective only; **no prospective P content was
  consumed or tested at PN4** (P's only prospective read remains
  H-GATE's single committed PN10 read, per P2 payload
  pattern_HGATE.prospective_test).
- **H-GATE status at this instance:** H-GATE (as stated in P2)
  concerns the ≤ 3 window; PN4's ≤ 3 window is EMPTY (re-anchored:
  457,450 supports, 0 nonzero — census_runner.w3_reanchor, plus the M
  mechanism theorem's corollary). The w=4 window is NONEMPTY at exactly
  d = 4 = min_i s_i, with 35 supports — matching Theorem T part 2
  (existence at weight d; NOTE: the predecessor called the first hit an
  instance of "part 1" — part 1 is the emptiness STRICTLY BELOW d;
  this is part 2) and the NEW exact characterization found here:
  at |S| = d the supports carrying codewords are EXACTLY the minimum
  lines of the argmin-d direction (all 35, no other supports of size 4
  carry any codeword). **Nothing falsifies H-GATE**; its prospective
  content (one committed read, PN10) stands confirmed in P2's payload
  and untouched by this campaign (we did not re-read PN10).
- **H-MINLINE note (owner's separate conjecture):** our 35=35 set
  equality at PN4 is CONSISTENT with it (its third pre-registered
  instance). Per the campaign contract it is NOT folded into W — status
  and evidence for H-MINLINE live in its own pre-registered campaign,
  not here. Recorded as "consistent; separate campaign owner".

## E. Defects disclosed (all in the additive instrument; none touched
frozen bytes, none touched the M campaign, none touched the LADDER)

1. `parity_crosscheck.py` w4 census, set-bug: the first
   `parity_w4_ledger.jsonl` crashed on JSON-serializing a Python set
   ("class" field). Process crash-looped 6 restarts; ledger burned
   [400k, 2.4M] of hit-before-progress bookkeeping. Preserved as
   `DEFECT_2026-09-01_parity_w4_ledger_set_bug.jsonl`. Fixed by sorting
   class vectors into lists before ledger write.
2. Same file, idx-increment bug (introduced by the fix edit): the
   recorder's `idx` never advanced → progress/hits/final records carried
   next_index 0 while the generator advanced (run continued through ~5M
   supports, harvest list correct; records misleading). Killed;
   preserved as `DEFECT2_2026-09-01_parity_w4_ledger_zero_counter.jsonl`;
   fixed by restoring `idx += 1`. The mathematical harvest in that run's
   memory was confirmed correct post-mortem (12 hits = direction-0 lines
   o=0..11 at their exact lex ids), so the formulation was never wrong —
   only the bookkeeping.
3. Both bugs were in MY additive w=4 code paths, not in the frozen
   pre-statement, not in census_runner.py (which we did not edit), not in
   the ladder. The final run (pid 87297, 04:46:57Z start) completed with
   zero defects: final next_index = 15,329,615, hits_total = 35, and the
   per-hit lex indices re-verify exactly.
4. During freeze prep: stray `__pycache__` (census_runner.pyc from the
   predecessor's run; parity_crosscheck.pyc from mine) — no fresh .pyc
   retained; both MOVED OUT of the campaign dir tree (preserved outside
   as checksummed evidence in DEFECT_NOTE_pyc_* names) and the empty
   dir removed, per the workspace control plane. census_runner.py was
   never modified (mtime unchanged, 1788231651).
5. The owner's instrument had a leaked-loop-variable bug earlier today
   (two false FAILs) — our treatment of owner INPUT as targets, never as
   evidence, was correct: all owner numbers used here (35 / set-equality
   / dims / d1/d2 line facts) were REPRODUCED by our own instruments
   before entering this record.

## F. Grandfathered legacy + rule-7 scope sentence (pre-registered domain)

**Scope (rule 7, same sentence shape as the pre-statement §5):** Swept:
the weight-EXACTLY-4 point-support window at (421,(4,5,7),(1,1,1),Λ=Id),
COMPLETE census of all 15,329,615 supports, PARITY-primary instrument
pre-registered, candidates = normalized reduced-basis of each nonzero
V∩F^S, per-candidate delta closed EXACTLY by the frozen exhaustive
ascending-cost line-tuple sweep, ρ^window_w4 = 1. NOT swept: every other
(q,s,t,Λ); weights ≠ 4 (the ≥5 and ≤3 windows keep their pre-existing
status); ρ_inst (still open over heavier supports; certified interval
now ρ_inst ∈ [1/420, 1/1] with 1 the w4-window cap... the ≤ 3-window
value is vacuous (window empty), the certified-instance interval is
ρ_inst ∈ [1/420, ρ^window_w4] = [1/420, 1] per the P2 convention, ρ_inst
itself NOT certified); nonprime fields; non-Id Λ (ρ is Λ-dependent per
the pre-statement's L7 note); H-MINLINE's other instances; the TR26-150
paper (untouched); H-GATE's PN10 prospective read (not re-run, status
inherited from P2's frozen payload).
**Grandfathered:** this campaign dir + pre-statement (commit 8114a8f)
predate the 2026-09-01 workspace campaign-lifecycle tool; per Main's
ruling it is finished in place, not migrated or re-registered.

## G. Artifacts (this freeze)

ADDITIVE files only (frozen pre_statement.md / census_runner.py /
anchors logs / kill-demo ledger untouched):
- `parity_crosscheck.py` (modified vs 11cc15a: w4-primary instrument +
  fixes; the as-run copy is what was executed by pid 87297) — sha256 in
  checksums_addendum.sha256
- `parity_w4_ledger.jsonl` (39 progress + 1 final; checksummed lines)
- `parity_w4_hits.json` (the 35 candidates, dims, classes, lex idx)
- `rho_w4.py`, `rho_w4_results.json` (35 exact delta closures; ρ^window_w4 = 1)
- `w4_vbasis_verification.json` (35/35 V-basis agreement on the found set)
- `w4_owner_input_checks.json` (our rerun of the owner-input line facts)
- `DEFECT_*`, `DEFECT2_*`, `DEFECT_NOTE_pyc_*` (disclosure evidence)
- NOT landed: B11/B12/B13 V-basis rows (process still running at freeze;
  see B) — census_ledger.jsonl intentionally NOT modified post-freeze.
- Named next action: **when w_ladder (pid 85270) drains B11-B13, append a
  one-line follow-up addendum with the V-basis full-space w=4 result at
  PN4 (or its honest non-landing) and the B11-B13 ladder tallies, and
  re-freeze with a checksums addendum2.** The ladder process itself is
  NOT killed by this freeze.
