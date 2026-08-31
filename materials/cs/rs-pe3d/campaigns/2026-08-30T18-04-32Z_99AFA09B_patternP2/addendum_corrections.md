# ADDENDUM — corrections recorded under rule 5 (inline retractions), 2026-08-30T23:5xZ

1. The manifest pattern_P_verdict string and the README originally said the
   PN weak arm landed 3/3. At freeze-build time PN6/PN7 were folded into
   the checkpoint but the verdict summary strings were written from the
   earlier landing set. Correct count: PN landed rows are PN1, PN2, PN3,
   PN6, PN7 = 5, all hits => 5/5 on landed weak-arm rows (same verdict,
   larger n). payload.json pattern_P.arm_PN rows_landed says 3 for the same
   stale reason: the authoritative per-row table (payload rows array) is
   correct (5 landing PN rows present); only the two summary strings were
   stale. Neither affects any computed value, witness, or the P verdict.

2. The README illustrative uniform-null figure was written as
   2.0e-7; the correct value of (1/5)^9 is 1/1953125 = 5.12e-7. The
   ILLUSTRATIVE label (not a rigorous null, not a p-value) stands.

3. This directory is append-only after the freeze; these corrections are
   recorded here unmodified, and payload.json/manifest.json/checksums are
   left as frozen. The README current-state section (a living document)
   is corrected in place with the same two fixes.
