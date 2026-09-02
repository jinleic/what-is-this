# Attempt log (freeze-hygiene: per-attempt distinct files)

- attempt 1 (attempt1_stdout/stderr.log): first two runs crashed at the
  rewritten comp-triple block (NameError cand — the min subtree got
  clobbered; then Ival-from-float in tier_calc). Fixed: restored the
  frozen parent min construction verbatim, built tier_calc arbs from
  exact arb endpoints of the recorded pairs. THIRD physical run:
  EXIT=0, all controls pass (E0, GATE-Z both records, C1, C2 region-2,
  C3 region-1 drift), P1/PV/P2 clean,
  VERDICT FROZEN-CERTIFIED, pinned: R_glob[2] (Tier-2 at BOTH e0 and e1:
  exclusion gap [4.5483061219230014e-12, 2.0511278540081266e-06], low
  4.55e-12 > 0 outward; the triplicate-identical e0==e1 endpoint values
  are float-render artifacts — the underlying exact intervals agree at
  render precision because the comp/glob candidate triples move by < 1
  in the last printed digit at r0=1e-7).
- Tier summary: R_comp[0,3]/[1,3]/[2,3]/[0,2] and R_glob[1] = Tier-3;
  R_glob[1]'s e1 exclusion gap lo = -2.54e-12 (straddles by ~2.5e-12 —
  a 62x-signal-lobe near-pin, recorded as the closest miss);
  R_glob[2] = Tier-2 EXCLUSION CERTIFIED (cand_2 excluded from the
  possible-minimizer set on the entire ray family; the argmin at the
  endpoint is {0,1} instead of {0,1,2}).
