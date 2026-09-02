# Route F run notes (interim — machine log digest, numbers to be re-derived
# from frozen artifacts before any verdict wording)

## Run 1 (broken, superseded): instrument ran the seeds against the RAW
(L_1, L_i, L_j) slices. The frozen seeds were optimized against the CONJUGATED
target tau ⊠ s (blockdiag doubled tau; both npz 'convention' fields state it,
and s3_routeF3_global.py V1 proves the T-similitude exactly). Symptom: seed
rel anchors read 8.17e-1 instead of 1.33e-3 / 1.00e-4. This is a target-
convention error in the instrument, not in the frozen data (rule 17c: a
correct measurement of the wrong object). Diagnostics from run 1 that are
TARGET-INDEPENDENT and stand:
  - 3-way table anchor: PASS (cd == gateA == upstream, unit, norm-mult).
  - rank(J) at the polished CP seed = 192/247 (precondition for Slice-1).
  - QRCP min pivot 1.75e-6 (float diagnostic >= 1e-9: Slice-1 fires).
  - exact rank(J_F) = 55/55 (full F-block transversality at that point).
  - 192x192 exact fmpq inverse: 34.5 s CPU single-threaded.

## Convention fix (verified before run 2)
Target rebuilt as blockdiag(tau_p, tau_p) with the entrywise similitude
anchor T L T == block (T = diag(1,...,1,-1)) asserted at startup; seed rel
anchors then reproduce EXACTLY:
  CP:  1.3252060344282049e-03, diff 0.00e+00 (bit-exact)
  EXT: 9.9987449458480189e-05, diff 6.16e-14 (float noise, within 1e-12
  trigger)
Both pass the pre-registered 1e-12 anchor. Run 2 (bg_14) is the certification
run of record.

## Known-honest expectation for the main system (rule 14 pre-thought)
The polish stalls near rel ~ 4e-5 after 5x5000 nfev (run 1 numbers on the
wrong target; run 2 numbers will differ but the stall pattern is structural:
LM/TRF averaged over 3 near-degenerate optical isomers). The Krawczyk window
needs ||Y||*||g0|| < rho << sigma_min(J_S)/(sum|factors|*132) — with
sigma_min ~ 1.75e-3-ish on the QRCP-selected block and |g0| ~ 4e-5 * ||T||,
the window likely does NOT open. Expected honest outcome: FAILURE TO CERTIFY
with certified exclusion/containment data per rung. Any PASS would have to
come from a rung where the polished residual is at the exact-arithmetic floor
(~1e-16 rel), which the declared 5-round budget does not reach. NO box
widening or seed additions are permitted to change that (Main condition 4;
rule 16).
