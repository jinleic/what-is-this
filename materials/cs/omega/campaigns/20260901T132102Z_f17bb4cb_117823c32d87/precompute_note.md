# Pre-compute notes (disclosed, append-only)

- attempted run 1 (attempt1_stdout/stderr.log): E0 PASS; crashed NameError T0
  (runner defect, engine untouched). No verdict.
- attempted run 2: crashed NameError SCRATCH (same class). No verdict.
- attempted run 3: OPEN at E-FIRST — exact r0 recomputed 0 vs expected
  1.0607e-6. ROOT CAUSE: runner's `inv_smax` comparison bug
  (`q < inv_smax` picks the MAX |s_i| inverse = 1.0, not the min) AND the
  positivity loop indexing bug (`c_frac[i]/(-s_exact[i])` with the tuple
  unpack shadow). Exact recomputation gives R_exact =
  4771454744492531/295147905179352825856 = 1.6166317499672092e-05 —
  the FROZEN CandidateWitness positivity_ceiling EXACTLY.
  REGISTRATION CORRECTION (pre-verdict, pre-P): the expected r0 numbers in
  pre_statement §2 were my own float transcription error — the FROZEN
  protocol_static_corrected.json positivity_ceiling is
  1.6166317499672092e-05 (not 1.06e-6, which was R/16 of the LP ray step,
- IMPORTANT correction to my own note line above ("1e-7 / max|s_i|"):
  the |s_i| box wall does NOT bind for this ray because the swept ray is
  certified against the positivity wall ONLY — s has max |s_i| = 38/32? No:
  re-derived exactly: max|s_i| over the 45 coordinates is <= 1 (s in
  {-1,0,+1}^45 per the frozen direction record), so the box wall is
  1e-7/1 = 1e-7 >> positivity 1.6166e-5? NO — 1.6166e-5 > 1e-7, so THE BOX
  WALL BINDS: R_exact = min(1e-7, 1.6166e-5) = 1e-7. Correction to the
  correction: the certified feasible radius of the ray inside BOTH the box
  (|delta_i| <= 1e-7) and the positivity cone is EXACTLY the registered
  box radius 1e-7 (ray endpoint = the frozen CandidateWitness step at
  divisor 1, radius [944473296573929, 9444732965739290427392] =
  1e-7 = 0x1.ad7f29abcaf48p-24 exactly as registered in §1). The
  positivity ceiling 1.6166e-5 is NOT binding at radius 1e-7.
  REGISTERED EXPECTATION (final, binding): r0 = 1e-7 exactly (the frozen
  binary64 radius), R_exact_ceiling_positivity = 1.6166317499672092e-05.
- attempt 2 (attempt2_stdout/stderr.log): E0/GATE-Z(±lemma1)/C1/C2/C3/P1/PV
  ALL PASS in-run (GATE-Z engine equivalence + identity accept + wrong-index
  reject + dual-drift reject + kernel membership of p*+r0*s all clean).
  CRASH in P2 leaf construction: frac_dyadic asserted q>0 on a ZERO delta
  coordinate — 7 of the 45 ray offsets s_i are 0, so abs(delta_i)=0 and a
  zero-width dyadic radius is the correct object. Fix: return the exact
  zero Ival for q == 0 (a point leaf at the center), before the q>0 assert.
  No verdict produced by attempt 2 (crash predates the P2 aggregation).
- attempt 3 (attempt3_stdout/stderr.log): all controls to P2 pass; P2 itself
  ran the full box aggregation; crash in the adjudication writer
  (report["P"] vs report["P2"] key typo, post-aggregation, engine clean).
  P2 aggregation was COMPLETED before the crash (9 s wall) — checkpoint
  saved through C1a. Fix is display-key only; no arithmetic touched.
- attempt 4 (attempt4_stdout/stderr.log): EXIT=0, VERDICT BOX-CONTRADICTION.
  All controls pass in-run: E0 (5 deps hash-verified), GATE-Z (with/without
  lemma1, frozen records reproduced), C1 (identity-accept, gap 5.79e-13),
  C2 (wrong-index reject, eps1 bit-exact vs frozen), C3 (dual-drift reject,
  coupling gap 4.41e-16 = 4*drift, branch response -2*prop*4*drift within
  2.9e-16), C1a (all non-dist0 leaves width-0; eps0 bit-constant at the ray
  endpoint 2.1502086942121845e-06), PV (kernel membership exact:
  A·delta==0, sum==0, |delta|<=r0=1e-7), PFIRST (walls exact: box binds).
- P (the pin, per Amendment-1 semantics): G1/G2 at e0 are the frozen point
  negative (−1.4335e-6 lows, the 2*prop_0*eps_0 charge floor) and at e1:
  G1 in [-1.8897456060889468e-05, +1.6396874550250076e-05] (width
  3.5294330611139544e-05 = 223.15x signal), G2 in
  [-1.9088262419827114e-05, +1.6587298910286668e-05] (width
  3.5675561330113786e-05 = 225.56x signal). Both straddle zero DECISIVELY:
  the pc-vs-pen order is NOT pinnable even on the exact-kernel ray family
  at full 1e-7 reach. Hull span cand0 at e1: 2.1109e-5 (vs parent's
  whole-box 3.2123e-5 — a 1.52x reduction only, still 133x the signal).
  VERDICT: BOX-CONTRADICTION (registered FROZEN-NEGATIVE branch).
- authoritative rerun in THIS dir after Amendment-2 commit 7305623:
  EXIT=0, decision fields identical (deterministic), 8.9 s wall.
